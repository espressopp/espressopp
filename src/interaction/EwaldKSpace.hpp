// ESPP_CLASS
#ifndef _INTERACTION_EWALDKSPACE_HPP
#define _INTERACTION_EWALDKSPACE_HPP

//#include <vector>
//#include <complex>

#include <boost/signals2.hpp>
#include "mpi.hpp"
#include "Potential.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"
#include "iterator/CellListIterator.hpp"
#include "Cell.hpp"
#include "bc/BC.hpp"

#include "System.hpp"

#include "boost/serialization/vector.hpp"
#include "boost/serialization/complex.hpp"

using namespace std;

typedef complex<double> dcomplex;

#define M_2PI (2*M_PIl)
#define M_PI2 (M_PIl*M_PIl)

#define M_1_SQRTPI (M_2_SQRTPIl * .5) /* 2/sqrt(pi)/2 = 1/sqrt(pi) */

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of the EwaldKSpace part*/
    class EwaldKSpace : public PotentialTemplate< EwaldKSpace > {
    private:
      real prefactor;
      real alpha; // Ewald parameter
      int kmax; // cutoff in k space
      // we need the system object to be able to access the box dimensions, communicator
      // number of particles, signals
      shared_ptr< System > system;
      
      real Lx, Ly, Lz; // system size the information could be got from bc
      int nParticles;  // number of particles 
      real rclx;
      real rcly;
      real rclz;
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

      // exponent array
      vector< vector<dcomplex> > eikx;
      vector< vector<dcomplex> > eiky;
      vector< vector<dcomplex> > eikz;
      
      vector< vector<dcomplex> > eik;
      
      real sum_q2;
        
      // @todo one should clean the memory then
      dcomplex* sum; 
      dcomplex* totsum; 
    public:
      static void registerPython();

      EwaldKSpace(shared_ptr< System > _system, real _prefactor, real _alpha, int _kmax);
      
      // at this point we are ready to prepare the kvector[], it can be done just once at the begin
      void preset(){
        // @todo it could be parallelized too !!!
        Real3D Li = system -> bc -> getBoxL(); // getting the system size
        Lx = Li[0];
        Ly = Li[1];
        Lz = Li[2];
        
        int kmaxsq = kmax * kmax;
        
        rclx = M_2PI / Lx;
        rcly = M_2PI / Ly;
        rclz = M_2PI / Lz;
        
        force_prefac[0] = prefactor * (-2.0) * M_2PI / Lx;
        force_prefac[1] = prefactor * (-2.0) * M_2PI / Ly;
        force_prefac[2] = prefactor * (-2.0) * M_2PI / Lz;
        
        // precalculate factors
        real B = M_PI2 / (alpha * alpha);
        real V = Lx * Ly * Lz * M_2PI;
        
        /* calculate the k-vector array */
        int ksq,kx,ky,kz;
        real rksq, rkx,rky,rkz;
        kVectorLength = 0;
        for(kx = 0; kx <= kmax; kx++){
          rkx = kx / Lx;
          for(ky = -kmax; ky <= kmax; ky++){
            rky = ky / Ly;
            for(kz = -kmax; kz <= kmax; kz++) {
              rkz = kz / Lz;
              ksq=kx*kx+ky*ky+kz*kz;
              if( (ksq < kmaxsq) && (ksq !=0) ){
                rksq = rkx*rkx+rky*rky+rkz*rkz;
                
                kvector.push_back(  exp( -rksq * B)/( rksq * V)  );
                
                kxfield.push_back( kx );
                kyfield.push_back( ky );
                kzfield.push_back( kz );
                
                kx_ind.push_back(        kx );
                ky_ind.push_back( kmax + ky );
                kz_ind.push_back( kmax + kz );
                
                kVectorLength++;
              }
            }
          }
        }
        
        sum = NULL;
        totsum = NULL;
        sum = new dcomplex[kVectorLength]; 
        totsum = new dcomplex[kVectorLength]; 
      }
      
      // here we get the current particle number on the current node
      // and set the auxiliary arrays eikx, eiky, eikz
      void getParticleNumber() {
        nParticles = system->storage->getNRealParticles();
        
        eikx = vector< vector<dcomplex> > (  kmax+1, vector<dcomplex>(nParticles, 0));
        eiky = vector< vector<dcomplex> > (2*kmax+1, vector<dcomplex>(nParticles, 0));
        eikz = vector< vector<dcomplex> > (2*kmax+1, vector<dcomplex>(nParticles, 0));
        
        eik = vector< vector<dcomplex> > (kVectorLength, vector<dcomplex>(nParticles, 0));
      }
      
      void count_charges(CellList realcells){
        real node_sum_q2 = 0.0;
        for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
          Particle &p = *it;
          node_sum_q2 += (p.q() * p.q());
        }
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
          sum[k] = dcomplex(0.0, 0.0);
          totsum[k] = dcomplex(0.0, 0.0); // @todo do we need to do that?
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
          node_energy += fact * kvector[k] *  norm( totsum[k] );
        }
        real energy = 0;
        mpi::all_reduce( *system -> comm, node_energy, energy, plus<real>() );
        
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


      // @todo this functions should be somehow cleaned
      real _computeEnergySqrRaw(real distSqr) const {
        cout << "This function currently doesn't work" << endl;
        return 0.0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        cout << "This function currently doesn't work" << endl;
        return false;
      }
      
    protected:
      boost::signals2::connection connectionRecalcKVec;
      boost::signals2::connection connectionGetParticleNumber;
    };
  }
}

#endif
