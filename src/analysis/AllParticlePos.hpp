// ESPP_CLASS
#ifndef _ANALYSIS_ALLPARTICLEPOS_HPP
#define _ANALYSIS_ALLPARTICLEPOS_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include <map>

namespace espresso {
  namespace analysis {

    /** Class provides pid and position information of all particles on all cpus
        (used e.g. for atom decomposition parallel analysis)
    */

    struct sBuf {
      real r[3];
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
   	    for (int i = 0; i < 3; ++i) ar & r[i];
       }
    };

    typedef std::map< size_t, sBuf > ConfMap;

    class AllParticlePos : public SystemAccess {
    public:
      AllParticlePos(shared_ptr<System> system) : SystemAccess (system) {};
      ~AllParticlePos() {};
      /** gather and broadcast all particle positions to all cpus */
      void gatherAllPositions();
      static void registerPython();

      ConfMap AllPositions;
      int numParticles;

    protected:
      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}

#endif
