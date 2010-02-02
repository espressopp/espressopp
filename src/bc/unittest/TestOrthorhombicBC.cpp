#define PARALLEL_TEST_MODULE BoundaryConditions
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "../OrthorhombicBC.hpp"

using namespace espresso;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr< bc::OrthorhombicBC > pbc;

  Fixture() {
    Real3D L(10.0, 10.0, 10.0);
    pbc = make_shared< bc::OrthorhombicBC >(L);
  }
};

BOOST_FIXTURE_TEST_CASE(foldingTest, Fixture) {
  int dim = 2;
  BOOST_CHECK_EQUAL(10.0, pbc->getBoxL(dim));

  Real3D pi(5.0, 5.0, 5.0);
  Real3D pj(11.0, 11.0, 11.0);
  Real3D rij;
  real d = 0.0;
  pbc->getMinimumImageVector(rij, d, pi, pj);
  BOOST_CHECK_EQUAL(rij[0], 4.0);
}
