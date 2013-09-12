#include "python.hpp"
#include "PressureTensorLayer.hpp"

namespace espresso {
  namespace analysis {

    void PressureTensorLayer::registerPython() {
      using namespace espresso::python;
      class_<PressureTensorLayer, bases< AnalysisBase > >
        ("analysis_PressureTensorLayer", init< shared_ptr< System >, real, real >())
        .add_property("h0",
              &PressureTensorLayer::getH0,
              &PressureTensorLayer::setH0)
        .add_property("dh",
              &PressureTensorLayer::getDH,
              &PressureTensorLayer::setDH)
      ;
    }
  }
}
