#include "python.hpp"
#include "PressureTensorMultiLayer.hpp"

namespace espresso {
  namespace analysis {

    void PressureTensorMultiLayer::registerPython() {
      using namespace espresso::python;
      class_<PressureTensorMultiLayer, bases< AnalysisBase > >
        ("analysis_PressureTensorMultiLayer", init< shared_ptr< System >, int, real >())
        .add_property("n",
              &PressureTensorMultiLayer::getN,
              &PressureTensorMultiLayer::setN)
        .add_property("dh",
              &PressureTensorMultiLayer::getDH,
              &PressureTensorMultiLayer::setDH)
      ;
    }
  }
}
