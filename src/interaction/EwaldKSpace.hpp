// ESPP_CLASS
#ifndef _INTERACTION_EWALDKSPACE_HPP
#define _INTERACTION_EWALDKSPACE_HPP

#include "Potential.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"
#include "iterator/CellListIterator.hpp"
#include "Cell.hpp"
#include "mpi.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of the EwaldKSpace part*/
    class EwaldKSpace : public PotentialTemplate< EwaldKSpace > {
    private:
      real alpha;
      int kmax;
      
      int * kvector;
    public:
      static void registerPython();

      // empty constructor the variables should be initialized later
      EwaldKSpace(): alpha(0.0), kmax(0) {
        preset();
      }
      EwaldKSpace(real _alpha, int _kmax): alpha(_alpha), kmax(_kmax) {
        preset();
      }

      // at this point we are ready to prepare the kvector[], it can be done just once at the begin
      void preset(){
//        /* measure the length of the array kvector[] */
//        int ksq,kx,ky,kz;
//        int totk=0;
//        for(kx = 0; kx <= ewald.kmax; kx++)
//          for(ky = -ewald.kmax; ky <= ewald.kmax; ky++) 
//            for(kz = -ewald.kmax; kz <= ewald.kmax; kz++) {
//              ksq=kx*kx+ky*ky+kz*kz;
//              if((ksq < ewald.kmaxsq) && (ksq !=0)) totk++;
//            }           
//          
//        double rkx,rky,rkz,rksq;
//        
//        
      }
      
      // set/get the parameters
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

      real _computeEnergy() const {
        real energy = 0;
        return energy;
      }

      bool _computeForce(CellList realcells) const {
        for (iterator::CellListIterator it(realcells); !it.isDone(); ++it) {
          Particle p = *it;
          //printf("Particle ID: %i\n",p.q());
          
          //mpi::all_reduce(mpi::, myforce, sumforce, boost::mpi::sum<real>());
        }
      }


      real _computeEnergySqrRaw(real distSqr) const {
        
      }

      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {

      }


    };
  }
}

#endif
