#ifndef _LBOUTPUT_PROFILE_VZOFX_HPP
#define _LBOUTPUT_PROFILE_VZOFX_HPP

#include "LBOutput.hpp"

namespace espresso {
  namespace analysis {
    class LBOutputProfileVzOfX : public LBOutput {
      public:
        LBOutputProfileVzOfX(shared_ptr<System> _system,
                          shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann);

        /** Destructor for output. */
/*        ~LBOutputProfileVzOfX ();
*/
        void writeOutput();

        static void registerPython();
    };
  }
}

#endif
