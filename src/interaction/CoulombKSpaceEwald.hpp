/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _INTERACTION_COULOMBKSPACEEWALD_HPP
#define _INTERACTION_COULOMBKSPACEEWALD_HPP

#include <cmath>

#include <boost/signals2.hpp>
#include "mpi.hpp"
#include "Potential.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"
#include "iterator/CellListIterator.hpp"
#include "Cell.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

#include "Tensor.hpp"

#include "System.hpp"

#include "boost/serialization/vector.hpp"
#include "boost/serialization/complex.hpp"

using namespace std;

typedef complex<espressopp::real> dcomplex;

// the following two constants are not defined everywhere (e.g. not in Mac OS X)
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_2_SQRTPIl
#define M_2_SQRTPIl 1.1283791670955125738961589031215452L
#endif

#define M_2PI (2*M_PIl)
#define M_PI2 (M_PIl*M_PIl)

#define M_1_SQRTPI (M_2_SQRTPIl * 0.5) /* 2/sqrt(pi)/2 = 1/sqrt(pi) */

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of the
     *  CoulombKSpaceEwald part. Currently it works with cubes and rectangular cuboids.
     *  Does not work for triclinic box, slab geometry.
     */
    
    // TODO should be optimized (force energy and virial calculate the same stuff)
    
    class CoulombKSpaceEwald : public PotentialTemplate< CoulombKSpaceEwald > {
    private:
      real prefactor;
      real alpha; // Ewald parameter
      int kmax; // cutoff in k space
      
      shared_ptr< System > system; // we need the system object to be able to access the box
                                   // dimensions, communicator, number of particles, signals
      
      real Lx, Ly, Lz; // local variable for system size
      int nParticles;  // local variable for the number of particles
      real rclx, rcly, rclz;
      real force_prefac[3]; // real array for force prefactors [0]: x,[1]: y,[2]: z
      
      vector<real> kvector; // precalculated k-vector
      int kVectorLength; // length of precalculated k-vector
      
      // This arrays and kVectorLength help to substitute 3 loops over kx, ky, kz by one
      vector<int> kxfield;
      vector<int> kyfield;
      vector<int> kzfield;
      // the same as kxfield, kyfield, kzfield, but adjusted as array indexes
      vector<int> kx_ind;
      vector<int> ky_ind;
      vector<int> kz_ind;
      
      // precalculated factors for the virial and virial tensor
      vector<real> virialPref;
      vector<Tensor> virialTensorPref;
      Tensor I;

      // exponent array
      vector< vector<dcomplex> > eikx;
      vector< vector<dcomplex> > eiky;
      vector< vector<dcomplex> > eikz;
      
      vector< vector<dcomplex> > eik;
      
      real sum_q2;
        
      dcomplex* sum;
      dcomplex* totsum;
    public:
      static void registerPython();

      CoulombKSpaceEwald(shared_ptr< System > _system, real _prefactor, real _alpha, int _kmax);
      
      ~CoulombKSpaceEwald();
      
      // at this point we are ready to prepare the kvector[], it can be done just once at the begin
      void preset(){
        // TODO it could be parallelized too
        Real3D Li = system -> bc -> getBoxL(); // getting the system size
        Lx = Li[0];
        Ly = Li[1];
        Lz = Li[2];
       
        real skmax = kmax / min(Lx, min(Ly,Lz)); 
        real skmaxsq = skmax * skmax ; // we choose the biggest cutoff 
        
        rclx = M_2PI / Lx;
        rcly = M_2PI / Ly;
        rclz = M_2PI / Lz;
        
        force_prefac[0] = prefactor * (-2.0) * rclx;
        force_prefac[1] = prefactor * (-2.0) * rcly;
        force_prefac[2] = prefactor * (-2.0) * rclz;
        
        // precalculate factors
        real invAlpha2 = 1.0 / (alpha * alpha);
        real B = M_PI2 * invAlpha2; // PI^2 / alpha^2
        real inv2alpha2 = 0.5 * invAlpha2; // 1.0 / (2*alpha^2)
        real V = M_2PI * Lx * Ly * Lz;
        
        /* calculate the k-vector array */
        int ksq, kx2, ky2, kz2;
        real rksq, rkx2, rky2, rkz2;
        real rLx2 = 1. / ( Lx * Lx );
        real rLy2 = 1. / ( Ly * Ly );
        real rLz2 = 1. / ( Lz * Lz );
        real rk2PIx, rk2PIy, rk2PIz;
        kVectorLength = 0;
        // clear all vectors
        kvector.clear();
        kxfield.clear();
        kyfield.clear();
        kzfield.clear();

        kx_ind.clear();
        ky_ind.clear();
        kz_ind.clear();
        virialPref.clear();
        virialTensorPref.clear();
        for(int kx = 0; kx <= kmax; kx++){
          kx2  = kx * kx;
          rkx2 = kx2 * rLx2;
          rk2PIx = kx * rclx;
          for(int ky = -kmax; ky <= kmax; ky++){
            ky2  = ky * ky;
            rky2 = ky2 * rLy2;
            rk2PIy = ky * rcly;
            for(int kz = -kmax; kz <= kmax; kz++) {
              kz2  = kz * kz;
              rkz2 = kz2 * rLz2;
              rk2PIz = kz * rclz;
              
              ksq  = kx2+ky2+kz2;
              rksq = rkx2+rky2+rkz2;
              
              if( (rksq < skmaxsq) && (ksq !=0) ){
                kvector.push_back(  exp( -rksq * B)/( rksq * V)  );
                
                kxfield.push_back( kx );
                kyfield.push_back( ky );
                kzfield.push_back( kz );
                
                kx_ind.push_back(        kx );
                ky_ind.push_back( kmax + ky );
                kz_ind.push_back( kmax + kz );
                
                // the tensor should be: deltaKronecker(i,j) - 2*hi*hj / h^2 - hi*hj / (2*alfa^2)
                Real3D h(rk2PIx, rk2PIy, rk2PIz);
                real h2 = h * h;
                Tensor hh(h, h); // it is tensor: hi*hj
                
                virialPref.push_back( 1 - h2 * inv2alpha2 );
                virialTensorPref.push_back( I - 2 * hh / h2 - hh * inv2alpha2 );
                
                kVectorLength++;
              }
            }
          }
        }
        
 //cout <<"node:  "<< system->comm->rank() <<  " kVectorLength: "<< kVectorLength<< "   kmax: "<< skmax  <<endl;
        if(sum != NULL){
          delete [] sum;
          sum = NULL;
        }
        if(totsum != NULL){
          delete [] totsum;
          totsum = NULL;
        }
        sum = new dcomplex[kVectorLength]; 
        totsum = new dcomplex[kVectorLength];
        
        getParticleNumber();
        
      }
      
      // here we get the current particle number on the current node
      // and set the auxiliary arrays eikx, eiky, eikz
      void getParticleNumber() {
        nParticles = system->storage->getNRealParticles();
        
        eikx = vector< vector<dcomplex> > (  kmax+1, vector<dcomplex>(nParticles, 0));
        eiky = vector< vector<dcomplex> > (2*kmax+1, vector<dcomplex>(nParticles, 0));
        eikz = vector< vector<dcomplex> > (2*kmax+1, vector<dcomplex>(nParticles, 0));
        
        eik  = vector< vector<dcomplex> > (kVectorLength, vector<dcomplex>(nParticles, 0));
      }
      
      // it counts the squared charges over all system. It is used for self energy calculations
      void count_charges(CellList realcells){
        real node_sum_q2 = 0.0;
        for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
          Particle &p = *it;
          node_sum_q2 += (p.q() * p.q());
        }
        sum_q2 = 0.0;
        mpi::all_reduce( *system -> comm, node_sum_q2, sum_q2, plus<real>() );
      }
      
      // set/get the parameters
      void setPrefactor(real _prefactor) {
      	prefactor = _prefactor;
        preset();
      }
      real getPrefactor() const { return prefactor; }
      void setAlpha(real _alpha) {
        alpha = _alpha;
        preset();
      }
      real getAlpha() const { return alpha; }
      void setKMax(int _kmax) {
        kmax = _kmax;
        preset();
      }
      int getKMax() const { return kmax; }

      // compute force and energy
      void exponentPrecalculation(CellList realcells){
        /* Calculation of k space sums */
        // -1, 0, 1
        int j=0; // auxiliary variable, particle counter
        for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
          Particle &p = *it;
          
          eikx[     0][j] = dcomplex(1.0, 0.0);
          eiky[kmax+0][j] = dcomplex(1.0, 0.0);
          eikz[kmax+0][j] = dcomplex(1.0, 0.0);

          eikx[     1][j] = dcomplex( cos( rclx * p.position()[0] ), sin( rclx * p.position()[0] ) );
          eiky[kmax+1][j] = dcomplex( cos( rcly * p.position()[1] ), sin( rcly * p.position()[1] ) );
          eikz[kmax+1][j] = dcomplex( cos( rclz * p.position()[2] ), sin( rclz * p.position()[2] ) );

          eiky[kmax-1][j] = conj( eiky[kmax+1][j] );
          eikz[kmax-1][j] = conj( eikz[kmax+1][j] );
         
          j++;
        }
        
        // calculation of the rest terms
        for (int k=2; k<=kmax; k++) {
          j=0;
          for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
            Particle &p = *it;

            eikx[     k][j] = eikx[     k-1][j] * eikx[     1][j];
            
            eiky[kmax+k][j] = eiky[kmax+k-1][j] * eiky[kmax+1][j];
            eiky[kmax-k][j] = conj( eiky[kmax+k][j] );

            eikz[kmax+k][j] = eikz[kmax+k-1][j] * eikz[kmax+1][j];
            eikz[kmax-k][j] = conj( eikz[kmax+k][j] );

            j++;
          }
        }
        
        real q; // particle charge, auxiliary
        for (int k=0; k<kVectorLength; k++) {
          sum[k]    = dcomplex(0.0, 0.0);
          totsum[k] = dcomplex(0.0, 0.0);
          j=0;
          for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
            Particle &p = *it;
            if( (q = p.q())!=0 ){
              eik[k][j] = q  *  ( eikx[ kx_ind[k] ][j] * eiky[ ky_ind[k] ][j] * eikz[ kz_ind[k] ][j] );
              sum[k] += eik[k][j];
            }
            j++;
          }
        }
        
        mpi::all_reduce( *system -> comm, sum, kVectorLength, totsum, plus<dcomplex>() );
      }
      
      real _computeEnergy(CellList realcells){
        // exponent array
        exponentPrecalculation(realcells);
        
        mpi::communicator communic = *system->comm;
        
        int n_nodes = communic.size();
        int this_node = communic.rank();
        
        // TODO it could be a problem if   n_nodes > kVectorLength
        int numk = kVectorLength / n_nodes + 1;
        int mink = this_node * numk;
        int maxk = mink + numk;
        if(maxk>kVectorLength) maxk = kVectorLength;
        real fact;
        real node_energy = 0;
        for(int k=mink; k<maxk; k++) {
          if (kxfield[k]==0) 
            fact=1.0;
          else
            fact=2.0;
          node_energy += fact * kvector[k] * norm( totsum[k] );
        } 
          //cout <<"node:  "<< this_node << "  node energy: "<< node_energy << "  fact: "<< fact<< " kmax: "<< kVectorLength <<endl;
//exit(0);
        real energy = 0;
        mpi::all_reduce( communic, node_energy, energy, plus<real>() );
        
        /* self energy correction */
        energy -= sum_q2 * alpha * M_1_SQRTPI;
        
        energy *= prefactor;
        
        return energy;
      }

      // @TODO this function could be void, 
      bool _computeForce(CellList realcells){
        // exponent array
        exponentPrecalculation(realcells);

        real fact; // factor due to the symmetry
        for (int k=0; k<kVectorLength; k++) {
          if ( kxfield[k] == 0)
            fact=1.0;
          else
            fact=2.0;
          
          dcomplex tff = fact * kvector[k] * totsum[k]; // auxiliary complex factor
          int j=0;
          for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
            Particle& p = *it;
            
            real tf = 0.0;
            if( p.q()!=0 ){
              tf= imag( tff  *  conj( eik[k][j] ) );
            
              p.force()[0] += force_prefac[0] * tf * kxfield[k];
              p.force()[1] += force_prefac[1] * tf * kyfield[k];
              p.force()[2] += force_prefac[2] * tf * kzfield[k];
            }
            j++;
          }
          
        }
        
        return true;
      }
      
      // compute virial for this interaction
      // (!note: all particle interaction contains only one potential)
      real _computeVirial(CellList realcells){
        // TODO it's exactly the same part as energy has
        // exponent array
        exponentPrecalculation(realcells);
        
        mpi::communicator communic = *system->comm;
        
        int n_nodes = communic.size();
        int this_node = communic.rank();
        
        // TODO it could be a problem if   n_nodes > kVectorLength
        int numk = kVectorLength / n_nodes + 1;
        int mink = this_node * numk;
        int maxk = mink + numk;
        if(maxk>kVectorLength) maxk = kVectorLength;
        real fact;
        real node_virial = 0;
        for(int k=mink; k<maxk; k++) {
          if (kxfield[k]==0) 
            fact=1.0;
          else
            fact=2.0;
          node_virial += fact * virialPref[k] * kvector[k] * norm( totsum[k] );
        }
        real virial = 0;
        mpi::all_reduce( communic, node_virial, virial, plus<real>() );
        
        return virial;
      }
      
      // compute virial Tensor for this interaction
      // (!note: all particle interaction contains only one potential)
      Tensor _computeVirialTensor(CellList realcells){
        // TODO it's exactly the same part as energy does
        // exponent array
        exponentPrecalculation(realcells);
        
        mpi::communicator communic = *system->comm;
        
        int n_nodes = communic.size();
        int this_node = communic.rank();
        
        // TODO it could be a problem if   n_nodes > kVectorLength
        int numk = kVectorLength / n_nodes + 1;
        int mink = this_node * numk;
        int maxk = mink + numk;
        if(maxk>kVectorLength) maxk = kVectorLength;
        real fact;
        Tensor node_virialTensor = 0;
        for(int k=mink; k<maxk; k++) {
          if (kxfield[k]==0) 
            fact=1.0;
          else
            fact=2.0;
          node_virialTensor += fact * kvector[k] * norm( totsum[k] ) * virialTensorPref[k];
        }
        
        Tensor virialTensor(0.0);
        mpi::all_reduce( communic, node_virialTensor, virialTensor, plus<Tensor>());

        // using boost::mpi::all_reduce or reduce is very slow for <Tensor>
        // as a suggestion, one could try the line below:
        // mpi::all_reduce( communic, (double*)&node_virialTensor, 6, (double*)&virialTensor, plus<double>());

        return virialTensor;
      }
      
      real _computeEnergySqrRaw(real distSqr) const {
        esutil::Error err(system->comm);
        stringstream msg;
        msg << "There is no sense to call this function for Ewald summation";
        err.setException( msg.str() );
        return 0.0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        esutil::Error err(system->comm);
        stringstream msg;
        msg << "There is no sense to call this function for Ewald summation";
        err.setException( msg.str() );
        return false;
      }
      
    protected:
      // it's responsible for the k vector recalculation when the box size changes
      boost::signals2::connection connectionRecalcKVec;
      // --||-- when the particle number is changed
      boost::signals2::connection connectionGetParticleNumber;
    };
  }
}

#endif
