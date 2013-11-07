#ifndef _LBOUTPUT_SCREEN_HPP
#define _LBOUTPUT_SCREEN_HPP

#include "LBOutput.hpp"

namespace espresso {
  namespace analysis {
    class LBOutputScreen : public LBOutput {
      public:
      LBOutputScreen(shared_ptr<System> _system,
                          shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann);

        /** Destructor for output. */
/*        ~LBOutputScreen ();
*/
        void writeOutput();

        static void registerPython();
    };
  }
}

#endif
