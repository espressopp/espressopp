#include "python.hpp"
#include "GammaVariate.hpp"
#include "esutil/RNG.hpp"

using namespace boost;

namespace espresso {
  namespace esutil {
    GammaVariate::GammaVariate(shared_ptr< RNG > _rng,
				 const int alpha,
				 const real beta)
      : Super(*(_rng->getBoostRNG()), DistType(alpha, beta)), rng(_rng)
    {}
    
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    GammaVariate::registerPython() {
      using namespace espresso::python;

      real (GammaVariate::*pyCall)()
      	= &GammaVariate::operator();

      class_< GammaVariate >("esutil_GammaVariate",
			      init< shared_ptr< RNG > >())
      	.def("__call__", pyCall)
      	;
    }
  }
}
