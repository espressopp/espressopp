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

typedef complex<double> dcomplex;

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

namespace espresso {
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
      
      // Brillouin zones for the optimal influence function
      static const int brillouin = 1;
      
      vector< vector<real> > precalc_interp_caf;
      
      vector< vector<real> > mesh_shift;
      
      vector< vector<real> > d_op; 
      
      vector< vector< vector<real> > >gf; // influence function
      
      //vector< vector< vector< dcomplex > > > QQQ;
      vector<  dcomplex > QQQ;
      
      // reference points in the lattice, needed for the charge assignment
      vector<Int3D> g_ca;
      // mesh contributions of the charged particles
      vector< vector< vector< vector< real > > > > q_l;
      
      int nParticles;  // number of particles in system
      Real3D sysL;     // system size
      real sumq_2, sum_q2; // squared sum of charges and sum of squared charges
      
      real af_coef[8][7][7]; // matrix of predefined assigned function coefficients
      
      
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
        precalc_interp_caf = vector< vector<real> > (P, vector<real>(2*interpolation+1, 0.0) );
        precalc_interpol_charge_assignment_f();
        
        mesh_shift = vector< vector<real> >(3, vector<real>() );
        d_op = vector< vector<real> >(3, vector<real>() );
        for(int i=0;i<3;i++){
          mesh_shift[i] = vector<real>(M[i], 0.0);
          d_op[i] = vector<real>(M[i], 0.0);
        }
        
        calc_m_shift();
        
        calc_differential_operator();
        
        gf = vector< vector<vector<real> > > (M[0], 
                     vector<vector<real> >(M[1], 
                     vector<real>(M[2], 0.0)) );
        
        calc_opt_influence_function();
        
        // -----------------------------------------
        // charge assignment
        /*
        QQQ = vector< vector<vector<dcomplex> > > (M[0], 
                      vector<vector<dcomplex> >(M[1], 
                             vector<dcomplex>(M[2], 0.0)) );
         */
        QQQ = vector<dcomplex>(M[0]*M[1]*M[2], 0.0);
        
        
        g_ca = vector<Int3D>(0, Int3D(0) );

        q_l = vector< vector<vector< vector<real> > > > (nParticles,
                      vector<vector< vector<real> > >(M[0], 
                             vector< vector<real> >(M[1], 
                                     vector<real>(M[2], 0.0))) );
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
          for(int j = 0; j<P; j++)
            precalc_interp_caf[j][i+interpolation] = asignment_f(x, j, P);
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
        
        real coef  = 2.0 * M[0]*M[1]*M[2] / (sysL[0]*sysL[1]);

        real denom;
        Real3D nom, D;
        Int3D N;
        for ( N[0] = 0; N[0] < M[0]; N[0]++){
          for ( N[1] = 0; N[1] < M[1]; N[1]++){
            for ( N[2] = 0; N[2] < M[2]; N[2]++){
              if ( N == Int3D(0) )
                gf[ N[0] ][ N[1] ][ N[2] ] = 0.0;
              else{
                aliasing_sum( N, &nom, &denom);
                for(int i=0; i<3; i++){
                  D[i] = d_op[i][N[i]]; 
                }

                real D2 = D.sqr();

                gf[ N[0] ][ N[1] ][ N[2] ] = 
                        (D2 > 1e-10) ? coef * ( D*nom ) / ( D2*denom*denom ) : 0.0;
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
      
      real aliasing_sum(Int3D N, Real3D *nominator, real *denominator){
        /*
         * taken from Deserno's code
         *  calculates the aliasing sums in the nominator and denominator of the 
           expression for the optimal influence function (see  Hockney/Eastwood: 
           8-22, p. 275).
           N: is the n-vector for which the aliasing sum is to be performed.
           *nominatorX,*nominatorY,*nominatorZ : x-, y-, and z-component of the 
                                                 aliasing sum.
           *denominator : aliasing sum in the denominator. */

        Real3D ms(mesh_shift[0][N[0]], mesh_shift[1][N[1]], mesh_shift[2][N[2]]);
        
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
        for(int i=0; i<3; i++) out[i] = dround(x[i]+0.5);
        return out;
      }
      
      
      real _computeEnergy(CellList realCells){
        real energy = 0;
        
        common_part(realCells, 1);
        
        int MM[3];
        for(int i=0; i<3; i++) MM[i] = M[i];
        
        fftw_complex *inE, *outE;
        fftw_plan pE;
        inE = (fftw_complex*) fftw_malloc( M[0]*M[1]*M[2] * sizeof(fftw_complex));
        outE = (fftw_complex*) fftw_malloc( M[0]*M[1]*M[2] * sizeof(fftw_complex));
        pE = fftw_plan_dft(3, MM, inE, outE, FFTW_FORWARD, FFTW_ESTIMATE);
        
        inE = reinterpret_cast<fftw_complex*>( &QQQ[0] );
        
        fftw_execute(pE);
        
        // stupid slow way of copying an array
        for(int i=0; i<M[0]; i++){
          for(int j=0; j<M[1]; j++){
            for(int k=0; k<M[2]; k++){
              int index = k + M[2] * (j + M[1] * i);
              QQQ[index] = dcomplex(outE[index][0], outE[index][1]);
            }
          }
        }
        
        //std::cout<<"segfalll."<<std::endl;
        //exit(0);
        
        // TODO fix bug
        fftw_destroy_plan(pE);
        inE = NULL; outE = NULL;
        fftw_free(inE); fftw_free(outE);
        
        energy = 0.0;
        for (int i=0; i<M[0]; i++){
          for (int j=0; j<M[1]; j++){
            for (int k=0; k<M[2]; k++){
              int index = k + M[2] * (j + M[1] * i);
              //energy += gf[i][j][k] * norm( QQQ[i][j][k] ) ;
              energy += gf[i][j][k] * norm( QQQ[index] ) ;
            }
          }
        }

        // TODO sysL[0]?? what about [1] and [2]?
        energy *= ( C_pref * sysL[0] / ((real)M[0]*(real)M[1]*(real)M[2]*4.0*M_PIl) );
        // self energy and net charge correction: 
        energy -= C_pref * ( sum_q2 * alpha / sqrt(M_PIl)  +  
                        sumq_2 * M_PIl / (2.0* sysL[0]*sysL[1]*sysL[2] *pow(alpha,2)) );
        
        return energy;
      }
      
      void common_part(CellList realCells, int iii){
        real _2interp = 2.0 * interpolation;
        
        //initialize();
        
        // TODO check variable assignshift and provide better logic
        int assignshift = M[0]- floor((real)(P-1)/2.0);
        
        real  modadd1, modadd2;
        // odd and even interpolation order
        switch (P) {
          case 2 : case 4 : case 6 : 
            { modadd1 = 0.5; modadd2 = -0.5;} break;
          case 1 :case 3 : case 5 : case 7 : 
            { modadd1 = 0.0; modadd2 =  0.5;} break;
        }

        // alternative reference to the arrays G[i][3] 
        Int3D Gi, arg;

        g_ca.clear();
        QQQ.clear();
        //QQQ.reserve(M[0]*M[1]*M[2]);
        QQQ = vector<dcomplex>(M[0]*M[1]*M[2], dcomplex(0.0));
        
        std::cout<< "QQQ.size(): "<< QQQ.size()<< "  i="<< iii << std::endl;
        
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
          g_ca.push_back( Gi );
          
          // Calculate the mesh based charges
          real T1,T2,T3;
          for (int j = 0; j < P; j++) {
            int xpos = (Gi[0] + j) % M[0];
            T1 = p.q() * precalc_interp_caf[j][arg[0]];
            for (int k = 0; k < P; k++) {
              int ypos = (Gi[1] + k) % M[1];
              T2 = T1 * precalc_interp_caf[k][arg[1]];
              for (int l = 0; l < P; l++) {
                int zpos = (Gi[2] + l) % M[2];
                T3 = T2 * precalc_interp_caf[l][arg[2]];

                // specific for force !!!!!!!!
                q_l[p.id()][xpos][ypos][zpos] = T3;
                
                int index = zpos + M[2] * (ypos + M[1] * xpos);
                
                QQQ[index] += dcomplex(T3, 0.0);
              }
            }
          }
        }
        
        /*
        if(iii==1)
          exit(0);
         */
         
        
      }

      // @TODO this function could be void, 
      bool _computeForce(CellList realCells){

        common_part(realCells, 0);
        
        int MM[3];
        for(int i=0; i<3; i++) MM[i] = M[i];
        fftw_complex *in, *out;
        fftw_plan p;
        in = (fftw_complex*) fftw_malloc( M[0]*M[1]*M[2] * sizeof(fftw_complex));
        out = (fftw_complex*) fftw_malloc( M[0]*M[1]*M[2] * sizeof(fftw_complex));
        p = fftw_plan_dft(3, MM, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        in = reinterpret_cast<fftw_complex*>( &QQQ[0] );
        
        fftw_execute(p); // repeat as needed 
        
        dcomplex *outdc = reinterpret_cast<dcomplex*>( out );
        QQQ.assign(outdc, outdc+M[0]*M[1]*M[2]);
        
        fftw_destroy_plan(p);
        in = NULL; out = NULL;
        fftw_free(in); fftw_free(out);
        
        /*
        vector<vector<vector<vector<dcomplex> > > > phi;
        phi = vector< vector< vector<vector<dcomplex> > > > (3, 
                      vector< vector<vector<dcomplex> > > (M[0], 
                              vector<vector<dcomplex> >(M[1], 
                                     vector<dcomplex>(M[2], dcomplex(0.0) ))));
         */
        vector<vector<dcomplex> > phi;
        phi = vector<vector<dcomplex> > (3, 
                     vector<dcomplex> (M[0]*M[1]*M[2], dcomplex(0.0) ));
        
        
        // Calculate the supporting arrays phi_?_??:
        for (int i=0; i<M[0]; i++){
          for (int j=0; j<M[1]; j++){
            for (int k=0; k<M[2]; k++) {  

              // new definition:
              int index = k + M[2] * (j + M[1] * i);
              for (int ii=0; ii<3; ii++){
                //phi[ii][i][j][k] =  d_op[ii][i] * gf[i][j][k] * 
                phi[ii][index] =  d_op[ii][i] * gf[i][j][k] * 
                        swap_complex( conj( QQQ[index] ) );
              }
            }
          }
        }
        
        

        fftw_plan pB;
        fftw_complex *inB, *outB;
        int MMM[3];
        for(int i=0;i<3;i++) MMM[i] = M[i];
        inB  = (fftw_complex*) fftw_malloc( M[0]*M[1]*M[2] * sizeof(fftw_complex));
        outB = (fftw_complex*) fftw_malloc( M[0]*M[1]*M[2] * sizeof(fftw_complex));
        pB = fftw_plan_dft(3, MMM, inB, outB, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for(int l=0; l<3; l++){
          /*
          // stupid slow way of copying an array
          for(int i=0; i<M[0]; i++){
            for(int j=0; j<M[1]; j++){
              for(int k=0; k<M[2]; k++){
                int index = k + M[2] * (j + M[1] * i);
                inB[index][0] = phi[l][i][j][k].real();
                inB[index][1] = phi[l][i][j][k].imag();
              }
            }
          }*/
          
          inB = reinterpret_cast<fftw_complex*>( &phi[l][0] );

          fftw_execute(pB);

          dcomplex *outdcB = reinterpret_cast<dcomplex*>( outB );
          phi[l].assign(outdcB, outdcB+M[0]*M[1]*M[2]);
          
          /*
          // stupid slow way of copying an array
          for(int i=0; i<M[0]; i++){
            for(int j=0; j<M[1]; j++){
              for(int k=0; k<M[2]; k++){
                int index = k + M[2] * (j + M[1] * i);
                phi[l][i][j][k] = dcomplex(outB[index][0], outB[index][1]);
              }
            }
          }*/
        }

        fftw_destroy_plan(pB);
        
        inB = NULL; outB = NULL;
        
        fftw_free(inB); fftw_free(outB);
        
        int iii = 0;
        for(iterator::CellListIterator it(realCells); it.isValid(); ++it){
          Particle &p = *it;
          
          Real3D fff(0.0);
          
          for (int i = 0; i < P; i++) {
            int xpos = (g_ca[iii][0] + i) % M[0];
            for (int j = 0; j < P; j++) {
              int ypos = (g_ca[iii][1] + j) % M[1];
              for (int k = 0; k < P; k++) {
                int zpos = (g_ca[iii][2] + k) % M[2];
                
                int index = zpos + M[2] * (ypos + M[1] * xpos);
                
                Real3D fff_add( phi[0][index].real(),
                                phi[1][index].real(),
                                phi[2][index].real());

                fff += C_pref * q_l[p.id()][xpos][ypos][zpos]  *  
                        fff_add / (real)(M[0]*M[1]*M[2]);
              }
            }
          }

          p.force() -= fff;
          iii++;
        }
        
        // usual return from espresso
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
      
      
      real asignment_f(real x, int k, int P) {
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
