#include "python.hpp"
#include "UniformOnSphere.hpp"
#include "RNG.hpp"

using namespace boost;

namespace espresso {
  namespace esutil {
    UniformOnSphere::UniformOnSphere(shared_ptr< RNG > _rng)
      : Super(*(_rng->getBoostRNG()), DistType(3)), rng(_rng)
    {}

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    UniformOnSphere::registerPython() {
      using namespace espresso::python;

      Real3D (UniformOnSphere::*pyCall)() 
      	= &UniformOnSphere::operator();

      class_< UniformOnSphere >("esutil_UniformOnSphere",
      				init< shared_ptr< RNG > >())
      	.def("__call__", pyCall)
      	;
    }
  }
}
