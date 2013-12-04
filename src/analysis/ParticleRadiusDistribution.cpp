#include "python.hpp"
#include "ParticleRadiusDistribution.hpp"

namespace espresso {
  namespace analysis {

    void ParticleRadiusDistribution::registerPython() {
      using namespace espresso::python;
      class_<ParticleRadiusDistribution, bases< AnalysisBase > >
        ("analysis_ParticleRadiusDistribution", init< shared_ptr< System > >())
      ;
    }
  }
}
