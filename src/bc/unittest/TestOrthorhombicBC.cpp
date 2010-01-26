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
    real L[3] = {10.0, 10.0, 10.0};
    pbc = make_shared< bc::OrthorhombicBC >(L);
  }
};

BOOST_FIXTURE_TEST_CASE(foldingTest, Fixture) {
  int dim = 2;
  BOOST_CHECK_EQUAL(10.0, pbc->getBoxL(dim));

  real pi[3] = {5.0, 5.0, 5.0};
  real pj[3] = {11.0, 11.0, 11.0};
  real rij[3];
  real d = 0.0;
  pbc->getMinimumImageVector(rij, d, pi, pj);
  BOOST_CHECK_EQUAL(rij[0], 4.0);
}
