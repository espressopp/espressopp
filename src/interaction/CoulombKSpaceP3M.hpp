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

//#include "Cell.hpp"

//#include "Tensor.hpp"

//#include "System.hpp"

//#include "boost/serialization/vector.hpp"
//#include "boost/serialization/complex.hpp"

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
     *  CoulombKSpaceP3M part. Currently it works with cubes and rectangular cuboids.
     *  Does not work for triclinic box and slab geometry.
     */
    
    // TODO should be optimized (force energy and virial calculate the same stuff)
    
    class CoulombKSpaceP3M : public PotentialTemplate< CoulombKSpaceP3M > {
    private:
      shared_ptr< System > system; // we need the system object to be able to access the box
                                   // dimensions, communicator, number of particles, signals
      
      real alpha; // Ewald splitting parameter
      Int3D M; // number of mesh-points
      int P; // charge assignment order
      real rc; // cutoff in real space
      real C_pref; // Coulomb prefactor
      real epsilon;
      
    public:
      static void registerPython();

      CoulombKSpaceP3M(shared_ptr< System > _system,
                       real _alpha,
                       Int3D _M,
                       int _P,
                       real _rcut,
                       real _coulomb_prefactor,
                       real _epsilon
              );
      
      ~CoulombKSpaceP3M();
      
      // 
      void preset(){

      }
      
      
      // set/get the parameters
      void setPrefactor(real _prefactor) {
      	C_pref = _prefactor;
        preset();
      }
      real getPrefactor() const { return C_pref; }
      
      
      // here we get the current particle number on the current node
      // and set the auxiliary arrays eikx, eiky, eikz
      void getParticleNumber() {
      }
      
      // it counts the squared charges over all system. It is used for self energy calculations
      void count_charges(CellList realcells){
      }
      
      // 
      
      real _computeEnergy(CellList realCells){
        real energy = 0;
        
        //mpi::all_reduce( communic, node_energy, energy, plus<real>() );
        /* self energy correction */
        //energy -= sum_q2 * alpha * M_1_SQRTPI;
        //energy *= prefactor;
        
        
        return energy;
      }

      // @TODO this function could be void, 
      bool _computeForce(CellList realCells){
        
        
        // usual return from espresso
        return true;
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
      
    protected:
      // it's responsible for the k vector recalculation when the box size changes
      boost::signals2::connection connectionRecalcKVec;
      // --||-- when the particle number is changed
      boost::signals2::connection connectionGetParticleNumber;
    };
  }
}

#endif
