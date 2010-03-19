#define BOOST_TEST_MODULE LennardJones
#include "ut.hpp"
#include "interaction/LennardJones.hpp"

using namespace espresso;
using namespace interaction;

BOOST_AUTO_TEST_CASE(LJNoCutoff)
{
  LennardJones lj(1.0, 1.0, infinity);
  // minimum
  BOOST_CHECK_CLOSE(lj.computeEnergy(pow(2.0, 1./6.)), -1.0, 0.01);

  // small at large values
  BOOST_CHECK_SMALL(lj.computeEnergy(2.5), 5.0);
  BOOST_CHECK_SMALL(lj.computeEnergy(10.0), 1.0);
}

BOOST_AUTO_TEST_CASE(LJ)
{
  LennardJones lj(1.0, 1.0, 2.0);

  // small at cutoff and shortly before
  BOOST_CHECK_SMALL(lj.computeEnergy(1.9999), 0.01);
  BOOST_CHECK_SMALL(lj.computeEnergy(2.0), 0.01);
  // zero after cutoff
  BOOST_CHECK_EQUAL(lj.computeEnergy(10.0), 0.0);
}
