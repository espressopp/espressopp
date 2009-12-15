#define PARALLEL_TEST_MODULE BoundaryConditions
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "../BC.hpp"

using namespace espresso;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  bc::BC::SelfPtr pbc;

  Fixture() {
    real L[3] = {10.0, 10.0, 10.0};
    pbc = make_shared< bc::BC >(L);
  }
};

BOOST_FIXTURE_TEST_CASE(foldingTest, Fixture) {
  BOOST_CHECK_EQUAL(1.0, 1.0);
}

