#define BOOST_TEST_MODULE PBC

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../PBC.hpp"

using namespace espresso;

struct Fixture {

    static const real size = 10.0;

    espresso::bc::PBC pbc;
    Fixture() : pbc(size) {}
    ~Fixture() {}
};

//Triplet stores the scalar values (si and sj) for the unfolded particle positions
//and the expected separation value (es). All three coordinates can be checked at
//once using the same values.
struct Triplet {
  real si;
  real sj;
  real es;
  void set(real _si, real _sj, real _es) {si=_si; sj=_sj; es=_es;}
};


/** Test 1: Try various pairs to see if correct rij is found */

BOOST_FIXTURE_TEST_CASE(Minimum_Image_Test, Fixture) {

  std::vector<Triplet> v;
  Triplet t;

  t.set(1.0, 9.0, 2.0);     v.push_back(t);
  t.set(9.0, 1.0, -2.0);    v.push_back(t);
  t.set(4.0, 5.0, -1.0);    v.push_back(t);
  t.set(-5.0, 9.0, -4.0);   v.push_back(t);
  t.set(-51.0, -9.0, -2.0); v.push_back(t);
  t.set(-34.0, 72.0, 4.0);  v.push_back(t);

  for(std::vector<Triplet>::iterator it = v.begin(); it != v.end(); it++)
    {
      //set the positions of particles i and j
      Real3D ri(it->si);
      Real3D rj(it->sj);

      //set the correct or expected minimum image distance
      Real3D rijExpected(it->es);

      //compute the minimum image distance using PBC::getDist
      Real3D rijComputed = pbc.getDist(ri, rj);

      //compare the expected and correct values
      BOOST_CHECK_CLOSE(rijExpected.getX(), rijComputed.getX(), 1.e-10f);
      BOOST_CHECK_CLOSE(rijExpected.getY(), rijComputed.getY(), 1.e-10f);
      BOOST_CHECK_CLOSE(rijExpected.getZ(), rijComputed.getZ(), 1.e-10f);
    }
}
