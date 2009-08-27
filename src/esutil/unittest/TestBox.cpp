#include <acconfig.hpp>
#define BOOST_TEST_MODULE box
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cmath>
#include "../Box.hpp"

using namespace espresso::esutil;

typedef Vector3D<float> Real3D;

BOOST_AUTO_TEST_CASE(box_construct) {
  {
    Box<float> box;
    BOOST_CHECK_EQUAL(box.getLeft(), Real3D(0.0));
    BOOST_CHECK_EQUAL(box.getRight(), Real3D(0.0));
    BOOST_CHECK_EQUAL(box.getExtend(), Real3D(0.0));
  }
  {
    Box<float> box(Real3D(42, 1, 100), Real3D(20, 30, 40));
    BOOST_CHECK_SMALL((box.getLeft() - Real3D(20, 1, 40)).abs(), 1e-4f);
    BOOST_CHECK_SMALL((box.getRight() - Real3D(42, 30, 100)).abs(), 1e-4f);
    BOOST_CHECK_SMALL((box.getExtend() - Real3D(22, 29, 60)).abs(), 1e-4f);
  }
}
