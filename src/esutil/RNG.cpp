#include "python.hpp"
#include "RNG.hpp"
#include "mpi.hpp"
#include "types.hpp"

using namespace boost;

namespace espresso {
  namespace esutil {

    RNG::RNG(long _seed): boostRNG(make_shared< RNGType >(_seed + mpiWorld->rank())),
            normalVariate(*boostRNG, normal_distribution< real >(0.0, 1.0)),
            uniformOnSphereVariate(*boostRNG, uniform_on_sphere< real, Real3D >(3))
    		//gammaVariate(*boostRNG, gamma_distribution< real >(1, 1.0)), //TODO this line is nonsense: alpha=1 is trivial
    {}

    void RNG::seed(long _seed) {
      // Seed the RNG for the given CPU
      boostRNG->seed(_seed + mpiWorld->rank());
    }

    real RNG::operator()() { 
      variate_generator< RNGType&, uniform_01<> > uni(*boostRNG, uniform_01<>());
      return uni();
    }

    int RNG::operator()(int N) { 
      uniform_smallint< int > uni_dist(0, N-1);
      variate_generator< RNGType&, uniform_smallint< int > > uni(*boostRNG, uni_dist);
      return uni();
    }

    real RNG::normal() {
      return normalVariate();
    }

    real RNG::gamma(unsigned int alpha) {
    	gamma_distribution< real > gamma_dist(alpha, 1.0); //scale parameter \beta=1.0
    	variate_generator< RNGType&, gamma_distribution< real > > gamma_var(*boostRNG, gamma_dist);
    	return gamma_var();
    }

    
    Real3D RNG::uniformOnSphere() {
      return uniformOnSphereVariate();
    }

    shared_ptr< RNGType > RNG::getBoostRNG() {
      return boostRNG;
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    RNG::registerPython() {
      using namespace espresso::python;

      real (RNG::*pyCall1)() = &RNG::operator();
      int (RNG::*pyCall2)(int) = &RNG::operator();


      class_< RNG >("esutil_RNG", init< boost::python::optional< long > >())
        .def("seed", &RNG::seed)
        .def("__call__", pyCall1)
        .def("__call__", pyCall2)
        .def("normal", &RNG::normal)
        .def("gamma", &RNG::gammaOf1)
        .def("gamma", &RNG::gamma)
        .def("uniformOnSphere", &RNG::uniformOnSphere);
    }
  }
}

