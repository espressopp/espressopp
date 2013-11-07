#include "python.hpp"
#include "LBOutput.hpp"

namespace espresso {
  namespace analysis {
//    LOG4ESPP_LOGGER(LBOutput::theLogger, "LBOutput");
    /* for abstract class there is no need to define anything here,
     * except of python registration */

    void LBOutput::registerPython() {
      using namespace espresso::python;

      class_<LBOutput, bases< AnalysisBase >, boost::noncopyable >
      ("analysis_LBOutput", no_init)

        .def("writeOutput", pure_virtual(&LBOutput::writeOutput))
      ;
    }
  }
}
