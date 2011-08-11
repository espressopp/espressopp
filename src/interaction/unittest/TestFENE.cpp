#define BOOST_TEST_MODULE FENE
#include "ut.hpp"
#include "interaction/FENE.hpp"
#include "Real3D.hpp"

using namespace espresso;
using namespace interaction;

BOOST_AUTO_TEST_CASE(FENE_noCutoff)
{
  FENE fene(30.0, 0.0, 1.5, infinity);

  BOOST_CHECK_CLOSE(fene.computeEnergy(0.97), 18.279, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(0.80), 11.296, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(1.40), 69.147, 0.01);
}

BOOST_AUTO_TEST_CASE(FENE_noArgsConstructor)
{
  /* 
  TODO: why  are next 4 lines failing?
  FENE fene();
  fene.setK(30.0);
  fene.setR0(0.0);
  fene.setRMax(1.5);
  */

  FENE fene(0.0, 0.0, 0.0, infinity);
  fene.setK(30.0);
  fene.setR0(0.0);
  fene.setRMax(1.5);

  BOOST_CHECK_CLOSE(fene.computeEnergy(0.97), 18.279, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(0.80), 11.296, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(1.40), 69.147, 0.01);
}

BOOST_AUTO_TEST_CASE(FENE_force)
{
  FENE fene(30.0, 0.0, 1.5, infinity);
  Real3D dist(0.9, 0.1, 0.1);
  Real3D f = fene.computeForce(dist);

  BOOST_CHECK_CLOSE(f[0], -42.782, 0.01);
  BOOST_CHECK_CLOSE(f[1], -4.7535, 0.01);
  BOOST_CHECK_CLOSE(f[2], -4.7535, 0.01);
}
