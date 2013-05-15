// ESPP_CLASS
#ifndef _ANALYSIS_INTRACHAINDISTSQ_HPP
#define _ANALYSIS_INTRACHAINDISTSQ_HPP

#include "python.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"
#include "FixedPairList.hpp"
#include "AllParticlePos.hpp"

namespace espresso {
  namespace analysis {

    /** calculate mean squared intra-chain-distances
    */

    class IntraChainDistSq : public AllParticlePos {
    public:
      shared_ptr<FixedPairList> fpl;
      IntraChainDistSq(shared_ptr<System> system, shared_ptr<FixedPairList> _fpl) : AllParticlePos(system), fpl(_fpl) {};
      ~IntraChainDistSq() {};
      python::list compute();
      static void registerPython();

    protected:
      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}

#endif
