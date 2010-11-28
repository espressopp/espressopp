#define BOOST_TEST_MODULE SoftCosine
#include "ut.hpp"
#include "interaction/SoftCosine.hpp"

using namespace espresso;
using namespace interaction;

BOOST_AUTO_TEST_CASE(SCNoCutoff)
{
  SoftCosine sc(1.0, 2.5);
  BOOST_CHECK_CLOSE(sc.computeEnergy(0.0), 2.0, 0.01);
  BOOST_CHECK_CLOSE(sc.computeEnergy(1.0), 1.3090169943, 0.01);
  BOOST_CHECK_CLOSE(sc.computeEnergy(2.5), 0.0, 0.01);
}

BOOST_AUTO_TEST_CASE(SC)
{
  SoftCosine sc(1.0, 2.5, 0.0);
  BOOST_CHECK_CLOSE(sc.computeEnergy(0.0), 2.0, 0.01);
  BOOST_CHECK_CLOSE(sc.computeEnergy(1.0), 1.3090169943, 0.01);
  //BOOST_CHECK_CLOSE(sc.computeEnergy(2.5), 0.0, 0.01);
}
