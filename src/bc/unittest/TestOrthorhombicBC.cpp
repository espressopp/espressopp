#define PARALLEL_TEST_MODULE BoundaryConditions
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "Real3D.hpp"
#include "esutil/RNG.hpp"
#include "bc/OrthorhombicBC.hpp"

using namespace espresso;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr< bc::BC > bc;

  Fixture() {
    Real3D L(10.0, 10.0, 10.0);
    shared_ptr< esutil::RNG > rng 
      = make_shared< esutil::RNG >();
    bc = make_shared< bc::OrthorhombicBC >(rng, L);
  }
};

BOOST_FIXTURE_TEST_CASE(foldingTest, Fixture) {
  BOOST_CHECK_EQUAL(Real3D(10.0), bc->getBoxL());

  Real3D pi(5.0, 5.0, 5.0);
  Real3D pj(11.0, 11.0, 11.0);
  Real3D rij;
  bc->getMinimumImageVector(rij, pi, pj);
  BOOST_CHECK_EQUAL(rij[0], 4.0);
}
