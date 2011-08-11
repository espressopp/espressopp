#include "python.hpp"
#include "NormalVariate.hpp"
#include "esutil/RNG.hpp"

using namespace boost;

namespace espresso {
  namespace esutil {
    NormalVariate::NormalVariate(shared_ptr< RNG > _rng,
				 const real mean, 
				 const real sigma) 
      : Super(*(_rng->getBoostRNG()), DistType(mean, sigma)), rng(_rng)
    {}
    
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    NormalVariate::registerPython() {
      using namespace espresso::python;

      real (NormalVariate::*pyCall)() 
      	= &NormalVariate::operator();

      class_< NormalVariate >("esutil_NormalVariate",
			      init< shared_ptr< RNG > >())
      	.def("__call__", pyCall)
      	;
    }
  }
}
