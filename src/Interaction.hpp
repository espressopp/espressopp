#ifndef _INTERACTION_HPP
#define _INTERACTION_HPP

#include "Particle.hpp"
#include "BC.hpp"
#include "logging.hpp"
#include "types.hpp"
#include "Storage.hpp"
#include "VerletList.hpp"

namespace espresso {
  namespace interaction {

    /** Interaction base class. */

    class Interaction {

    public:

      // full loop over a storage; this class provides its own implementation
      // but will call the abstract routine for two Cells.

      virtual real computeStorageEnergy(shared_ptr<Storage> storage);

      virtual void addVerletListForces(shared_ptr<VerletList> vl) = 0;

      virtual real computeVerletListEnergy(shared_ptr<VerletList> vl) = 0;

      virtual real computeCellEnergy(ParticleList &pl) = 0;

      virtual real computeCellEnergy(ParticleList &pl1, ParticleList &pl2) = 0;

    public:

      class ParametersBase {
      private:
        real cutoff;
        real cutoffSqr;

      public:
        ParametersBase() { setCutoff(0.0); }
        void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff * cutoff; }
        real getCutoff() const { return cutoff; }
        real getCutoffSqr() const { return cutoffSqr; }
      };
    
    protected:

      /** Logger */

      static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#define INTERACTION_ROUTINES(Derived) \
 \
void Derived::addVerletListForces(shared_ptr<VerletList> vl) \
{ \
  addVerletListForcesImpl(vl); \
} \
 \
real Derived::computeVerletListEnergy(shared_ptr<VerletList> vl) \
{ \
  return computeVerletListEnergyImpl(vl); \
} \
 \
real Derived::computeCellEnergy(ParticleList &pl) \
{ \
  return computeCellEnergyImpl(pl); \
} \
 \
real Derived::computeCellEnergy(ParticleList &pl1, ParticleList &pl2) \
{ \
  return computeCellEnergyImpl(pl1, pl2); \
} 


#endif
