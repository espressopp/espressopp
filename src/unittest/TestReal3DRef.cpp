#define BOOST_TEST_MODULE Real3DRef

#include "ut.hpp"
#include "Real3D.hpp"

using namespace espresso;

BOOST_AUTO_TEST_CASE(FromCArray) {
  real v[3];
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;

  ConstReal3DRef cr(v);
  // check that the Real3DRef reflects the given values
  BOOST_CHECK_EQUAL(cr[0], 1.0);
  BOOST_CHECK_EQUAL(cr[1], 2.0);
  BOOST_CHECK_EQUAL(cr[2], 3.0);

  Real3DRef r(v);
  // check that the Real3DRef reflects the given values
  BOOST_CHECK_EQUAL(r[0], 1.0);
  BOOST_CHECK_EQUAL(r[1], 2.0);
  BOOST_CHECK_EQUAL(r[2], 3.0);

  // check that the Real3DRef can be written to
  r[0]=42.0;
  r[1]=52.0;
  r[2]=62.0;

  // check that wrinting to the Real3DRef modifies the original array
  BOOST_CHECK_EQUAL(v[0], 42.0);
  BOOST_CHECK_EQUAL(v[1], 52.0);
  BOOST_CHECK_EQUAL(v[2], 62.0);

}

BOOST_AUTO_TEST_CASE(FromReal3D) {
  Real3D v;
  
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;

  ConstReal3DRef cr(v);
  // check that the Real3DRef reflects the given values
  BOOST_CHECK_EQUAL(cr[0], 1.0);
  BOOST_CHECK_EQUAL(cr[1], 2.0);
  BOOST_CHECK_EQUAL(cr[2], 3.0);

  Real3DRef r(v);
  // check that the Real3DRef reflects the given values
  BOOST_CHECK_EQUAL(r[0], 1.0);
  BOOST_CHECK_EQUAL(r[1], 2.0);
  BOOST_CHECK_EQUAL(r[2], 3.0);
 

  // check that the Real3DRef can be written to
  r[0]=42.0;
  r[1]=52.0;
  r[2]=62.0;

  // check that writing to the Real3DRef modifies the original
  BOOST_CHECK_EQUAL(v[0], 42.0);
  BOOST_CHECK_EQUAL(v[1], 52.0);
  BOOST_CHECK_EQUAL(v[2], 62.0);

}

BOOST_AUTO_TEST_CASE(AssignmentIsNotInit) {
  Real3D v(1.0, 1.0, 1.0);
  Real3D v2(2.0, 2.0, 2.0);

  Real3DRef r(v);
  BOOST_CHECK_EQUAL(r[0], 1.0);
  BOOST_CHECK_EQUAL(r[1], 1.0);
  BOOST_CHECK_EQUAL(r[2], 1.0);
  
  r = v2;
  BOOST_CHECK_EQUAL(r[0], 2.0);
  BOOST_CHECK_EQUAL(r[1], 2.0);
  BOOST_CHECK_EQUAL(r[2], 2.0);

  r[0] = 3.0;
  r[1] = 3.0;
  r[2] = 3.0;

  // Now, v should have changed, but v2 shouldn't!
  BOOST_CHECK_EQUAL(v[0], 3.0);
  BOOST_CHECK_EQUAL(v[1], 3.0);
  BOOST_CHECK_EQUAL(v[2], 3.0);
  BOOST_CHECK_EQUAL(v2[0], 2.0);
  BOOST_CHECK_EQUAL(v2[1], 2.0);
  BOOST_CHECK_EQUAL(v2[2], 2.0);
}

BOOST_AUTO_TEST_CASE(at) {
  Real3D a(1.0, 2.0, 3.0);

  ConstReal3DRef cr(a);
  BOOST_CHECK_EQUAL(cr.at(0), 1.0);
  BOOST_CHECK_EQUAL(cr.at(1), 2.0);
  BOOST_CHECK_EQUAL(cr.at(2), 3.0);
  BOOST_CHECK_THROW(cr.at(-1), std::out_of_range);
  BOOST_CHECK_THROW(cr.at(3), std::out_of_range);

  Real3DRef r(a);
  BOOST_CHECK_EQUAL(r.at(0), 1.0);
  BOOST_CHECK_EQUAL(r.at(1), 2.0);
  BOOST_CHECK_EQUAL(r.at(2), 3.0);
  BOOST_CHECK_THROW(r.at(-1), std::out_of_range);
  BOOST_CHECK_THROW(r.at(3), std::out_of_range);

  r.at(0) = 42.0;
  r.at(1) = 52.0;
  r.at(2) = 62.0;

  BOOST_CHECK_EQUAL(r[0], 42.0);
  BOOST_CHECK_EQUAL(r[1], 52.0);
  BOOST_CHECK_EQUAL(r[2], 62.0);
}

BOOST_AUTO_TEST_CASE(constElementAccess) {
  Real3D a(1.0, 2.0, 3.0);
  const Real3DRef r(a);

  BOOST_CHECK_EQUAL(r[0], 1.0);
  BOOST_CHECK_EQUAL(r[1], 2.0);
  BOOST_CHECK_EQUAL(r[2], 3.0);

  BOOST_CHECK_EQUAL(r.at(0), 1.0);
  BOOST_CHECK_EQUAL(r.at(1), 2.0);
  BOOST_CHECK_EQUAL(r.at(2), 3.0);
  BOOST_CHECK_THROW(r.at(-1), std::out_of_range);
  BOOST_CHECK_THROW(r.at(3), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(testSqr) {
  Real3D a(1.0, 2.0, 3.0);
  Real3DRef r(a);
  BOOST_CHECK_CLOSE(r.sqr(), 14.0, 0.0001);

  ConstReal3DRef cr(a);
  BOOST_CHECK_CLOSE(cr.sqr(), 14.0, 0.0001);
}
  
BOOST_AUTO_TEST_CASE(testAbs) {
  Real3D a(1.0, 2.0, 3.0);
  Real3DRef r(a);
  BOOST_CHECK_CLOSE(r.abs(), sqrt(14.0), 0.0001);

  ConstReal3DRef cr(a);
  BOOST_CHECK_CLOSE(cr.abs(), sqrt(14.0), 0.0001);
}

BOOST_AUTO_TEST_CASE(unaryOps) {
  Real3D a(1.0, 2.0, 3.0);
  Real3D b(42.0, 52.0, 62.0);
  Real3DRef r(a);
  Real3DRef s(b);

  r += s;
  BOOST_CHECK_CLOSE(r[0], 43.0, 0.0001);
  BOOST_CHECK_CLOSE(r[1], 54.0, 0.0001);
  BOOST_CHECK_CLOSE(r[2], 65.0, 0.0001);

  r -= s;
  BOOST_CHECK_CLOSE(r[0], 1.0, 0.0001);
  BOOST_CHECK_CLOSE(r[1], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(r[2], 3.0, 0.0001);

  r *= 2.0;
  BOOST_CHECK_CLOSE(r[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(r[1], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(r[2], 6.0, 0.0001);
  
  r /= 2.0;
  BOOST_CHECK_CLOSE(r[0], 1.0, 0.0001);
  BOOST_CHECK_CLOSE(r[1], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(r[2], 3.0, 0.0001);

  ConstReal3DRef cs(b);
  r += cs;
  BOOST_CHECK_CLOSE(r[0], 43.0, 0.0001);
  BOOST_CHECK_CLOSE(r[1], 54.0, 0.0001);
  BOOST_CHECK_CLOSE(r[2], 65.0, 0.0001);

  r -= cs;
  BOOST_CHECK_CLOSE(r[0], 1.0, 0.0001);
  BOOST_CHECK_CLOSE(r[1], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(r[2], 3.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(boolOps) {
  Real3D a(1.0, 2.0, 3.0);
  Real3D b(42.0, 52.0, 62.0);
  Real3D c(1.0, 2.0, 3.0);

  Real3DRef r(a);
  Real3DRef s(b);
  Real3DRef t(c);

  BOOST_CHECK_EQUAL(r, r);
  BOOST_CHECK(!(r!=r));

  BOOST_CHECK_NE(r, s);
  BOOST_CHECK(!(r==s));

  BOOST_CHECK_EQUAL(r, t);
  BOOST_CHECK(!(r!=t));

  ConstReal3DRef cr(a);
  ConstReal3DRef cs(b);
  ConstReal3DRef ct(c);
  
  BOOST_CHECK_EQUAL(cr, cr);
  BOOST_CHECK(!(cr!=cr));

  BOOST_CHECK_NE(cr, cs);
  BOOST_CHECK(!(cr==cs));

  BOOST_CHECK_EQUAL(cr, ct);
  BOOST_CHECK(!(cr!=ct));

  // cross comparison
  BOOST_CHECK_EQUAL(cr, r);
  BOOST_CHECK(!(cr!=r));
}

