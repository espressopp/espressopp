#define BOOST_TEST_MODULE Morse
#include "ut.hpp"
#include "interaction/Morse.hpp"

using namespace espresso;
using namespace interaction;

BOOST_AUTO_TEST_CASE(MRNoCutoff)
{
  Morse mr(1.0, 1.0, 1.0, infinity);
  // minimum is U(r) at r=rMin
  BOOST_CHECK_CLOSE(mr.computeEnergy(1.0), -1.0, 0.01);

  // large value
  BOOST_CHECK_SMALL(mr.computeEnergy(10.0), 0.001);
}

BOOST_AUTO_TEST_CASE(MR)
{
  Morse mr(1.0, 1.0, 1.0, 2.0);

  // small at cutoff and shortly before
  BOOST_CHECK_SMALL(mr.computeEnergy(1.9999), 0.01);
  BOOST_CHECK_SMALL(mr.computeEnergy(2.0), 0.01);
  // zero after cutoff
  BOOST_CHECK_EQUAL(mr.computeEnergy(2.0001), 0.0);
}
