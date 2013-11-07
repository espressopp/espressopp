// ESPP_CLASS
#ifndef _ANALYSIS_LBOUTPUT_HPP
#define _ANALYSIS_LBOUTPUT_HPP

#include "integrator/LatticeBoltzmann.hpp"
#include "types.hpp"
#include "AnalysisBase.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {
  namespace analysis {
    using namespace iterator;
    /** Abstract base class for arbitrary output from LB simulations. */
    class LBOutput : public AnalysisBaseTemplate< real > {
    public:
      /* Constructor for the class */
      LBOutput(shared_ptr< System > _system,
          shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann)
      : AnalysisBaseTemplate< real >(_system) {
        latticeboltzmann = _latticeboltzmann;
      }
      /* Destructor for the class */
      virtual ~LBOutput () {}

      real computeRaw() {
        writeOutput();
        real _dummy = 0.;
        return _dummy;
      }

      python::list compute() {
        python::list _dummy;
        real res = computeRaw();
        _dummy.append(0.);
        return _dummy;
      }

      python::list getAverageValue() {
        python::list _dummy;
        _dummy.append(0.);
        return _dummy;
      }

      void resetAverage() {}
      void updateAverage(real res) {return;}

      /* writing of a profile into the output */
      virtual void writeOutput () = 0;

      static void registerPython();
    protected:
      shared_ptr<integrator::LatticeBoltzmann> latticeboltzmann;
    };
  }
}

#endif
