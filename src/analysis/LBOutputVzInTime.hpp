#ifndef _LBOUTPUT_VZINTIME_HPP
#define _LBOUTPUT_VZINTIME_HPP

#include "LBOutput.hpp"

namespace espresso {
  namespace analysis {
    class LBOutputVzInTime : public LBOutput {
      public:
      LBOutputVzInTime(shared_ptr<System> _system,
                          shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann);

        /** Destructor for output. */
/*        ~LBOutputVzInTime ();
*/
        void writeOutput();

        static void registerPython();
    };
  }
}

#endif
