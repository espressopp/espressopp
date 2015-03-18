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
#ifndef _INTERACTION_COULOMBKSPACEP3M_HPP
#define _INTERACTION_COULOMBKSPACEP3M_HPP

#include <cmath>
#include <boost/signals2.hpp>

#include <fftw3.h>

#include "mpi.hpp"
#include "Potential.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"

#include "bc/BC.hpp"

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
     *  CoulombKSpaceP3M part.
     * 
     *  The code is based on M.Deserno's work. Reference in literature
     *  M. Deserno, C.Holm, J.Chem. Phys, 109[18] (1998) 7694
     */
    
    // TODO this version is still under development. not parallel. doesn't work for several
    // CPUs
    // TODO should be optimized (force, energy and virial calculate the same stuff)
    
    class CoulombKSpaceP3M : public PotentialTemplate< CoulombKSpaceP3M > {
    private:
      shared_ptr< System > system; // we need the system object to be able to access the box
                                   // dimensions, communicator, number of particles, signals
      
      Real3D d_mesh; // distance between two meshpoints
      Real3D d_inv; // inverted distance between two meshpoints
      
      real C_pref; // Coulomb prefactor
      real alpha; // Ewald splitting parameter
      Int3D M; // number of mesh-points
      int P; // charge assignment order
      real rc; // cutoff in real space
      int interpolation; // number of interpolation points for the charge assignment 
                        // function
      
      int MMM;  // MMM = M[0]*M[1]*M[2]
      
      // Brillouin zones for the optimal influence function
      static const int brillouin = 1;
      
      vector< vector<real> > precalc_interp_caf;
      
      vector< vector<real> > mesh_shift;
      
      vector< vector<real> > d_op; 
      
      //vector< vector< vector<real> > >gf; // influence function
      vector<real> gf; // influence function
      
      //vector< vector< vector< dcomplex > > > QQQ;
      vector<  dcomplex > QQQ;
      
      vector<vector<dcomplex> > phi;
        
      // reference points in the lattice, needed for the charge assignment
      vector<Int3D> g_ca;
      // mesh contributions of the charged particles
      vector< vector< real > > q_l;
      
      vector< vector< vector< int > > > map_indx;
      
      int nParticles;  // number of particles in system
      Real3D sysL;     // system size
      real sumq_2, sum_q2; // squared sum of charges and sum of squared charges
      
      real af_coef[8][7][7]; // matrix of predefined assigned function coefficients
      
      // fftw elements
      fftw_complex *in_array;
      fftw_plan plan;
      int MM[3]; // auxiliary array for plan generation
      
      //real oddeven1, oddeven2; // supporting variables odd/even interpolation order
    public:
      static void registerPython();

      CoulombKSpaceP3M(shared_ptr< System > _system,
                       real _coulomb_prefactor,
                       real _alpha,
                       Int3D _M,
                       int _P,
                       real _rcut,
                       int _interpolation
              );
      
      ~CoulombKSpaceP3M();
      
      // 
      void preset(){
        sysL = system -> bc -> getBoxL();
        MMM = M[0] * M[1] * M[2];
        for(int i=0; i<3; i++) MM[i] = M[i];
        
        // map xpos, ypos and zpos to index = zpos + M[2] * (ypos + M[1] * xpos)
        map_indx = vector< vector< vector<int> > > (M[0],
                           vector< vector<int> >   (M[1],
                                   vector<int>     (M[2], 0) ));
        for(int i=0; i<M[0]; i++){
          int mi = M[1] * i;
          for(int j=0; j<M[1]; j++){
            int mj = M[2] * (j + mi);
            for(int k=0; k<M[2]; k++){
              map_indx[i][j][k] = k + mj;
            }
          }
        }
        
        precalc_interp_caf = vector< vector<real> > (P, vector<real>(2*interpolation+1, 0.0) );
        precalc_interpol_charge_assignment_f();
        
      }
      
/////////////////////////////////////////////////////////////////////////////////////////
      //  setters and getters
      void setPrefactor(real _prefactor) {
      	C_pref = _prefactor;
        preset();
      }
      real getPrefactor() const { return C_pref; }
      void setAlpha(real _alpha) {
      	alpha = _alpha;
        preset();
      }
      real getAlpha() const { return alpha; }
      void setMesh(Int3D _M) {
      	M = _M;
        preset();
      }
      Int3D getMesh() const { return M; }
      void setP(int _P) {
      	P = _P;
        preset();
      }
      int getP() const { return P; }
      void setCutoff(real _rc) {
      	rc = _rc;
        preset();
      }
      real getCutoff() const { return rc; }
      void setInterpolation(int _interpolation) {
      	interpolation = _interpolation;
        preset();
      }
      int getInterpolation() const { return interpolation; }
/////////////////////////////////////////////////////////////////////////////////////////

      void initialize(){
        
        set_fftw_array();
        
        mesh_shift = vector< vector<real> >(3, vector<real>() );
        d_op = vector< vector<real> >(3, vector<real>() );
        for(int i=0;i<3;i++){
          mesh_shift[i] = vector<real>(M[i], 0.0);
          d_op[i] = vector<real>(M[i], 0.0);
        }
        
        calc_m_shift();
        
        calc_differential_operator();
        
        gf = vector<real>(MMM, 0.0);
        
        calc_opt_influence_function();
        
        // -----------------------------------------
        // charge assignment
        QQQ = vector<dcomplex>(MMM, 0.0);
        
        
        // TODO check double use in common part
        g_ca = vector<Int3D>(0, Int3D(0) );

        q_l = vector< vector<real> > (nParticles, vector<real>(MMM, 0.0) );
        
        // force specific
        phi = vector<vector<dcomplex> > (3, vector<dcomplex> (MMM, dcomplex(0.0) ));
      }
      
      // get the current particle number on the current node
      // and set the auxiliary arrays
      void getParticleNumber() {
        nParticles = system->storage->getNRealParticles();
      }
      
      // it counts the squared charges over all system. It is used for self energy calculations
      void count_charges(CellList realCells){
        real node_sumq_2, node_sum_q2;
        node_sumq_2 = node_sum_q2 = 0.0;
        for(iterator::CellListIterator it(realCells); it.isValid(); ++it){
          Particle &p = *it;
          node_sumq_2  += p.q();
          node_sum_q2  += pow( p.q(), 2 );
        }
        
        sumq_2 = sum_q2 = 0.0;
        mpi::all_reduce( *system -> comm, node_sumq_2, sumq_2, plus<real>() );
        mpi::all_reduce( *system -> comm, node_sum_q2, sum_q2, plus<real>() );
        
        sumq_2 *= sumq_2;
      }
      
      // 
      void precalc_interpol_charge_assignment_f(){
        real _2interpol_inv = 1.0 / (2.0*interpolation);
        
        for (int i=-interpolation; i<=interpolation; i++) {
          real x = i * _2interpol_inv;
          for(int j = 0; j<P; j++){
            //precalc_interp_caf[j][i+interpolation] = asignment_f2(j, x, P);
            precalc_interp_caf[j][i+interpolation] = asignment_f1(x, j, P);
          }
        }
      }
      
      // calculate differential operator
      void calc_differential_operator(){
        for(int i=0; i<3; i++){
          for(int j = 1; j < M[i]/2; j++) {
            d_op[i][j] = j;
            d_op[i][M[i] - j] = -j;
          }
          d_op[i][M[i]/2] = 0.0;
        }
      }

      void create_mesh(CellList realCells){
      }
      void gen_mesh(CellList realCells){
      }
      void assign_charge_for_single_particle(real q, Real3D particle_pos){
      }
      
      // calculates the optimal influence function
      void calc_opt_influence_function(){
        
        real coef  = 2.0 * MMM / (sysL[0]*sysL[1]);

        real denom;
        Real3D nom, D;
        Int3D i;
        for ( i[0] = 0; i[0] < M[0]; i[0]++){
          for ( i[1] = 0; i[1] < M[1]; i[1]++){
            for ( i[2] = 0; i[2] < M[2]; i[2]++){
              int indx = map_indx[i[0]][i[1]][i[2]];
              if ( i == Int3D(0) )
                gf[ indx ] = 0.0;
              else{
                aliasing_sum( i, &nom, &denom);
                for(int l=0; l<3; l++) D[l] = d_op[l][i[l]];
                real D2 = D.sqr();
                gf[ indx ] = (D2 > 1e-10) ? coef * ( D*nom ) / ( D2*denom*denom ) : 0.0;
              }
            }
          }
        }
        
      }
      
      void calc_m_shift(){
        for(int i=0; i<3; i++){
          for(int j = 1; j < M[i]/2; j++) {
            mesh_shift[i][j] = j;
            mesh_shift[i][M[i] - j] = -j;
          }
          mesh_shift[i][M[i]/2] = 0.0;
        }
      }
      
      void aliasing_sum(Int3D i, Real3D *nominator, real *denominator){
        /*
         * taken from Deserno's code
         *  calculates the aliasing sums in the nominator and denominator of the 
           expression for the optimal influence function (see  Hockney/Eastwood: 
           8-22, p. 275).
           N: is the n-vector for which the aliasing sum is to be performed.
           *nominatorX,*nominatorY,*nominatorZ : x-, y-, and z-component of the 
                                                 aliasing sum.
           *denominator : aliasing sum in the denominator. */

        Real3D ms(mesh_shift[0][i[0]], mesh_shift[1][i[1]], mesh_shift[2][i[2]]);
        
        real sc;
        Real3D nml(0.0);
        
        Real3D Minv;
        for(int i=0;i<3;i++) Minv[i] = 1.0/(real)M[i];
        
        real ftr = pow( M_PIl / (alpha*sysL[0]), 2);

        *denominator = 0.0;
        *nominator = Real3D(0.0);
        
        Int3D ijk(0);
        Real3D x(0.0);
        for ( ijk[0] = -brillouin; ijk[0] <= brillouin; ijk[0]++) {
          for ( ijk[1] = -brillouin; ijk[1] <= brillouin; ijk[1]++) {
            for ( ijk[2] = -brillouin; ijk[2] <= brillouin; ijk[2]++) {
              
              for(int i=0; i<3; i++){
                nml[i] = ms[i] + M[i] * ijk[i];
                x[i]   = Minv[i] * nml[i];                
              }

              sc  = pow( sinc(x), 2.0*P);

              *denominator += sc;

              real nml2 = nml.sqr();
              real e = ftr * nml2;
              
              *nominator += ( ( e<30 ) ? sc*exp(-e)/nml2 : 0.0 ) * nml;
            }
          }
        }
      }
      

      /**
       *  taken from Deserno's code
       *  Calculates the sinc-function as sin(PI*x)/(PI*x).
       *
       * (same convention as in Hockney/Eastwood). In order to avoid
       * divisions by 0, arguments, whose modulus is smaller than epsi, will
       * be evaluated by an 8th order Taylor expansion of the sinc
       * function. Note that the difference between sinc(x) and this
       * expansion is smaller than 0.235e-12, if x is smaller than 0.1. (The
       * next term in the expansion is the 10th order contribution
       * PI^10/39916800 * x^10 = 0.2346...*x^12).  This expansion should
       * also save time, since it reduces the number of function calls to
       * sin().  
      */
      real sinc(real x){
        real epsi = 0.1;

        real c2 = -0.1666666666667e-0;
        real c4 =  0.8333333333333e-2;
        real c6 = -0.1984126984127e-3;
        real c8 =  0.2755731922399e-5;

        real PIx = M_PIl*x, PIx2;

        if (fabs(x)>epsi){
          return sin(PIx)/PIx;
        }
        else {
          PIx2 = PIx * PIx;
          return 1.0 + PIx2*(c2+PIx2*(c4+PIx2*(c6+PIx2*c8)));
        }
      }
      real sinc(Real3D X){
        real res = 1.0;
        for(int i=0;i<3;i++) res *= sinc(X[i]);
        return res;
      }
      
      real dround(real x) { return floor(x+0.5); }
      Real3D dround(Real3D x) {
        Real3D out;
        for(int i=0; i<3; i++) out[i] = dround(x[i]); //+0.5
        return out;
      }
      
      void set_fftw_array(){
        in_array = (fftw_complex*) fftw_malloc( MMM * sizeof(fftw_complex));
      }
      void set_plan_frw(){
        plan = fftw_plan_dft(3, MM, in_array, in_array, FFTW_FORWARD, FFTW_ESTIMATE);
      }
      void set_plan_bcw(){
        plan = fftw_plan_dft(3, MM, in_array, in_array, FFTW_BACKWARD, FFTW_ESTIMATE);
      }
      void clean_fftw(){
        fftw_destroy_plan(plan);
        in_array = NULL;
        fftw_free(in_array);
      }
      
      
      real _computeEnergy(CellList realCells){
        
        common_part(realCells, 1);
        
        real energy = 0.0;
        
        for (int i=0; i<MMM; i++){
          energy += gf[i] * norm( QQQ[i] );
        }
        
        // TODO sysL[0]?? what about [1] and [2]?
        energy *= ( C_pref * sysL[0] / (4.0*MMM*M_PIl) );

        // self energy and net charge correction: 
        energy -= C_pref * ( sum_q2 * alpha / sqrt(M_PIl)  +  
                        sumq_2 * M_PIl / (2.0* sysL[0]*sysL[1]*sysL[2] *pow(alpha,2)) );
        
        return energy;
      }
      
      // TODO get rid of iii at the end
      void common_part(CellList realCells, int iii){
        initialize();
        
        real _2interp = 2.0 * interpolation;
        // TODO assignshift probably should be [3]
        int assignshift = M[0]- floor((real)(P-1)/2.0);
        
        real  modadd1, modadd2;
        // odd and even interpolation order
        switch (P) {
          case 2 : case 4 : case 6 : 
            { modadd1 = 0.5; modadd2 = -0.5;} break;
          case 1 :case 3 : case 5 : case 7 : 
            { modadd1 = 0.0; modadd2 =  0.5;} break;
        }

        g_ca.clear();
        g_ca = vector<Int3D>(nParticles, Int3D(0) );
        QQQ.clear();
        QQQ = vector<dcomplex>(MMM, dcomplex(0.0));
        
        Int3D Gi, arg;
        for(iterator::CellListIterator it(realCells); it.isValid(); ++it){
          Particle &p = *it;
          Real3D ppos = p.position();
          
          Real3D d1;
          for(int i=0; i<3; i++){
            d1[i] = ppos[i] * M[i] / sysL[i] + modadd1;
          }
          Gi  = Int3D(d1 + modadd2) + assignshift;
          arg = Int3D( (d1 - dround(d1) + 0.5)*_2interp );

          // specific for force !!!!!!!!
          g_ca[p.id()] = Gi;
          
          // Calculate the mesh based charges
          real T1,T2,T3;
          for (int i = 0; i < P; i++) {
            int xpos = (Gi[0] + i) % M[0];
            T1 = p.q() * precalc_interp_caf[i][arg[0]];
            for (int j = 0; j < P; j++) {
              int ypos = (Gi[1] + j) % M[1];
              T2 = T1 * precalc_interp_caf[j][arg[1]];
              for (int k = 0; k < P; k++) {
                int zpos = (Gi[2] + k) % M[2];
                T3 = T2 * precalc_interp_caf[k][arg[2]];
                
                int indx = map_indx[xpos][ypos][zpos];
                
                // specific for force !!!!!!!!
                q_l[p.id()][indx] = T3;
                
                QQQ[indx] += dcomplex(T3, 0.0);
              }
            }
          }

        }
 
        in_array = reinterpret_cast<fftw_complex*>( &QQQ[0] );
        set_plan_frw();
        fftw_execute(plan);
      }

      // @TODO this function could be void, 
      bool _computeForce(CellList realCells){

        common_part(realCells, 0);
        
        // Calculate the supporting arrays phi_?_??:
        Int3D i;
        for ( i[0]=0; i[0]<M[0]; i[0]++){
          for ( i[1]=0; i[1]<M[1]; i[1]++){
            for ( i[2]=0; i[2]<M[2]; i[2]++) {  
              int indx = map_indx[i[0]][i[1]][i[2]];
              dcomplex phi_aux = gf[indx] * swap_complex( conj( QQQ[indx] ) );
              
              for (int ii=0; ii<3; ii++) phi[ii][indx] = d_op[ii][i[ii]] * phi_aux;
            }
          }
        }

        for(int l=0; l<3; l++){
          in_array = reinterpret_cast<fftw_complex*>( &phi[l][0] );
          set_plan_bcw();
          fftw_execute(plan);
        }
        
        real C_MMM_inv = C_pref / (real)MMM;
        for(iterator::CellListIterator it(realCells); it.isValid(); ++it){
          Particle &p = *it;
          
          int iii = p.id();
          Real3D ff(0.0);
          for (int i = 0; i < P; i++) {
            int xpos = (g_ca[iii][0] + i) % M[0];
            for (int j = 0; j < P; j++) {
              int ypos = (g_ca[iii][1] + j) % M[1];
              for (int k = 0; k < P; k++) {
                int zpos = (g_ca[iii][2] + k) % M[2];
                
                int indx = map_indx[xpos][ypos][zpos];
                
                Real3D f_add( phi[0][indx].real(),
                              phi[1][indx].real(),
                              phi[2][indx].real());

                ff += C_MMM_inv * q_l[iii][indx]  *  f_add ;
              }
            }
          }

          p.force() -= ff;
        }
        
        // usual return from espressopp
        return true;
      }
    
      dcomplex swap_complex(dcomplex C){
        return dcomplex(C.imag(), C.real());
      }
      

      // compute virial for this interaction
      // (!note: all particle interaction contains only one potential)
      real _computeVirial(CellList realcells){
        real virial = 0;
        //mpi::all_reduce( communic, node_virial, virial, plus<real>() );
        
        return virial;
      }
      
      // compute virial Tensor for this interaction
      // (!note: all particle interaction contains only one potential)
      Tensor _computeVirialTensor(CellList realcells){
        
        Tensor virialTensor(0.0);
        //mpi::all_reduce( communic, node_virialTensor, virialTensor, plus<Tensor>());

        // using boost::mpi::all_reduce or reduce is very slow for <Tensor>
        // as a suggestion, one could try the line below:
        // mpi::all_reduce( communic, (double*)&node_virialTensor, 6, (double*)&virialTensor, plus<double>());

        return virialTensor;
      }
      
      real _computeEnergySqrRaw(real distSqr) const {
        esutil::Error err(system->comm);
        stringstream msg;
        msg << "There is no sense to call this function for P3M";
        err.setException( msg.str() );
        return 0.0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        esutil::Error err(system->comm);
        stringstream msg;
        msg << "There is no sense to call this function for P3M";
        err.setException( msg.str() );
        return false;
      }
      
      
      real asignment_f1(real x, int k, int P) {
        esutil::Error err(system->comm);
        stringstream msg;
        if(P<1 || P>7){
          msg<<"Charge assignment order P="<<P<<"  For some reason order "<<k<<" is used";
          err.setException( msg.str() );
        }
        err.checkException();
        
        real res = af_coef[P][k][0];
        real xx = 1.0;
        for(int i=1; i<P; i++){
          xx *= x;
          res += af_coef[P][k][i] * xx;
        }
        return res;
      }
       
      real asignment_f2(int i, real x, int P) {
      mpi::communicator communic = *system->comm;
      int this_node = communic.rank();
      switch (P) {
        case 1 : return 1.0;
        case 2 : {
          switch (i) {
          case 0: return 0.5-x;
          case 1: return 0.5+x;
          default:
            fprintf(stderr,"%d: Tried to access charge assignment "
                    "function of degree %d in scheme of order %d.\n",this_node,i,P);
            return 0.0;
          }
        }
        case 3 : {
          switch (i) {
          case 0: return 0.5 * (0.5 - x);
          case 1: return 0.75 - pow(x, 2);
          case 2: return 0.5 * pow(0.5 + x, 2);
          default:
            fprintf(stderr,"%d: Tried to access charge assignment function"
                    " of degree %d in scheme of order %d.\n",this_node,i,P);
            return 0.0;
          }
        case 4 : {
          switch (i) {
          case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
          case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
          case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
          case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
          default:
            fprintf(stderr,"%d: Tried to access charge assignment function"
                    " of degree %d in scheme of order %d.\n",this_node,i,P);
            return 0.0;
          }
        }
        case 5 : {
          switch (i) {
          case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
          case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
          case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
          case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
          case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
          default:
            fprintf(stderr,"%d: Tried to access charge assignment function"
                    " of degree %d in scheme of order %d.\n",this_node,i,P);
            return 0.0;
          }
        }
        case 6 : {
          switch (i) {
          case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
          case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
          case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
          case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
          case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
          case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
          default:
            fprintf(stderr,"%d: Tried to access charge assignment function"
                    " of degree %d in scheme of order %d.\n", this_node, i, P);
            return 0.0;
          }
        }
        case 7 : {
          switch (i) {
          case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
          case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
          case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
          case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
          case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
          case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
          case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
          default:
            fprintf(stderr,"%d: Tried to access charge assignment function"
                    " of degree %d in scheme of order %d.\n", this_node, i, P);
            return 0.0;
          }
        }
        default :{
          fprintf(stderr,"%d: Charge assignment order %d unknown.\n",this_node, P);
          return 0.0;
        }
        }
      }
    }
                                                                  
      
    protected:
      // it's responsible for the k vector recalculation when the box size changes
      boost::signals2::connection connectionRecalcKVec;
      // --||-- when the particle number is changed
      boost::signals2::connection connectionGetParticleNumber;
      
      //  ???before force calculation???
      //boost::signals2::connection recalcCommonPart;
    };
  }
}

#endif
