#include <acconfig.hpp>
#define BOOST_TEST_MODULE box
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cmath>
#include "../Box.hpp"

using namespace espresso::esutil;

typedef SmallVector<float, 2> Real2D;
typedef Vector3D<float> Real3D;

BOOST_AUTO_TEST_CASE(box_construct_2d) {
  {
    Box<Real2D> box;
    BOOST_CHECK_EQUAL(box.getLeft(), Real2D(0.0));
    BOOST_CHECK_EQUAL(box.getRight(), Real2D(0.0));
    BOOST_CHECK_EQUAL(box.getExtend(), Real2D(0.0));
  }
  {
    Real2D cornerIL = (float []){42,  1};
    Real2D cornerIR = (float []){20, 30};
    Real2D cornerOL = (float []){20,  1};
    Real2D cornerOR = (float []){42, 30};
    Real2D extend   = (float []){22, 29};

    Box<Real2D> box(cornerIL, cornerIR);

    BOOST_CHECK_SMALL((box.getLeft() - cornerOL).abs(), 1e-4f);
    BOOST_CHECK_SMALL((box.getRight() - cornerOR).abs(), 1e-4f);
    BOOST_CHECK_SMALL((box.getExtend() - extend).abs(), 1e-4f);
  }
}

BOOST_AUTO_TEST_CASE(box_construct_3d) {
  {
    Box<Real3D> box;
    BOOST_CHECK_EQUAL(box.getLeft(), Real3D(0.0));
    BOOST_CHECK_EQUAL(box.getRight(), Real3D(0.0));
    BOOST_CHECK_EQUAL(box.getExtend(), Real3D(0.0));
  }
  {
    Real3D cornerIL(42,  1, 100);
    Real3D cornerIR(20, 30, 100);
    Real3D cornerOL(20,  1, 100);
    Real3D cornerOR(42, 30, 100);
    Real3D extend  (22, 29,   0);

    Box<Real3D> box(cornerIL, cornerIR);

    BOOST_CHECK_SMALL((box.getLeft() - cornerOL).abs(), 1e-4f);
    BOOST_CHECK_SMALL((box.getRight() - cornerOR).abs(), 1e-4f);
    BOOST_CHECK_SMALL((box.getExtend() - extend).abs(), 1e-4f);
  }
}
