#ifndef _INTERACTION_HPP
#define _INTERACTION_HPP

#include "types.hpp"
#include "logging.hpp"
#include "esutil/ESPPIterator.hpp"

namespace espresso {
  namespace interaction {

    /** Interaction base class. */

    class Interaction {

    public:
      virtual void addForces() = 0;
      virtual real computeEnergy() = 0;

      // move into separate classes:
      // virtual void addStorageForces(shared_ptr< Storage > storage) = 0;
      // virtual real computeStorageEnergy(shared_ptr< Storage > storage) = 0;

      // virtual void addVerletListForces(shared_ptr< VerletList > vl) = 0;
      // virtual real computeVerletListEnergy(shared_ptr< VerletList > vl) = 0;

      /** This method returns the maximal cutoff defined for one type pair. */
      virtual real getMaxCutoff() = 0;

    protected:
      /** Logger */
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    struct InteractionList
      : public std::vector< shared_ptr< Interaction > > {
      typedef esutil::ESPPIterator< std::vector< Interaction > > Iterator;
    };

    /** Interaction function base class. */
    class InteractionFunction {
    private:
      real cutoff;
      real cutoffSqr;
      
    public:
      InteractionFunction() { setCutoff(0.0); }
      void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff * cutoff; }
      real getCutoff() const { return cutoff; }
      real getCutoffSqr() const { return cutoffSqr; }
    };

  }
}

#endif
