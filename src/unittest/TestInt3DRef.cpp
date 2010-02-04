#define BOOST_TEST_MODULE Int3DRef

#include "ut.hpp"
#include "Int3D.hpp"

using namespace espresso;

BOOST_AUTO_TEST_CASE(FromCArray) {
  int v[3];
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;

  ConstInt3DRef cr(v);
  // check that the Int3DRef reflects the given values
  BOOST_CHECK_EQUAL(cr[0], 1);
  BOOST_CHECK_EQUAL(cr[1], 2);
  BOOST_CHECK_EQUAL(cr[2], 3);

  Int3DRef r(v);
  // check that the Int3DRef reflects the given values
  BOOST_CHECK_EQUAL(r[0], 1);
  BOOST_CHECK_EQUAL(r[1], 2);
  BOOST_CHECK_EQUAL(r[2], 3);

  // check that the Int3DRef can be written to
  r[0]=42;
  r[1]=52;
  r[2]=62;

  // check that wrinting to the Int3DRef modifies the original array
  BOOST_CHECK_EQUAL(v[0], 42);
  BOOST_CHECK_EQUAL(v[1], 52);
  BOOST_CHECK_EQUAL(v[2], 62);
}

BOOST_AUTO_TEST_CASE(FromInt3D) {
  Int3D v;
  
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;

  ConstInt3DRef cr(v);
  // check that the Int3DRef reflects the given values
  BOOST_CHECK_EQUAL(cr[0], 1);
  BOOST_CHECK_EQUAL(cr[1], 2);
  BOOST_CHECK_EQUAL(cr[2], 3);

  Int3DRef r(v);
  // check that the Int3DRef reflects the given values
  BOOST_CHECK_EQUAL(r[0], 1);
  BOOST_CHECK_EQUAL(r[1], 2);
  BOOST_CHECK_EQUAL(r[2], 3);
 
  // check that the Int3DRef can be written to
  r[0]=42;
  r[1]=52;
  r[2]=62;

  // check that writing to the Int3DRef modifies the original
  BOOST_CHECK_EQUAL(v[0], 42);
  BOOST_CHECK_EQUAL(v[1], 52);
  BOOST_CHECK_EQUAL(v[2], 62);
}

BOOST_AUTO_TEST_CASE(AssignmentIsNotInit) {
  Int3D v(1, 1, 1);
  Int3D v2(2, 2, 2);

  Int3DRef r(v);
  BOOST_CHECK_EQUAL(r[0], 1);
  BOOST_CHECK_EQUAL(r[1], 1);
  BOOST_CHECK_EQUAL(r[2], 1);
  
  r = v2;
  BOOST_CHECK_EQUAL(r[0], 2);
  BOOST_CHECK_EQUAL(r[1], 2);
  BOOST_CHECK_EQUAL(r[2], 2);

  r[0] = 3;
  r[1] = 3;
  r[2] = 3;

  // Now, v should have changed, but v2 shouldn't!
  BOOST_CHECK_EQUAL(v[0], 3);
  BOOST_CHECK_EQUAL(v[1], 3);
  BOOST_CHECK_EQUAL(v[2], 3);
  BOOST_CHECK_EQUAL(v2[0], 2);
  BOOST_CHECK_EQUAL(v2[1], 2);
  BOOST_CHECK_EQUAL(v2[2], 2);
}

BOOST_AUTO_TEST_CASE(at) {
  Int3D a(1, 2, 3);

  ConstInt3DRef cr(a);
  BOOST_CHECK_EQUAL(cr.at(0), 1);
  BOOST_CHECK_EQUAL(cr.at(1), 2);
  BOOST_CHECK_EQUAL(cr.at(2), 3);
  BOOST_CHECK_THROW(cr.at(-1), std::out_of_range);
  BOOST_CHECK_THROW(cr.at(3), std::out_of_range);

  Int3DRef r(a);
  BOOST_CHECK_EQUAL(r.at(0), 1);
  BOOST_CHECK_EQUAL(r.at(1), 2);
  BOOST_CHECK_EQUAL(r.at(2), 3);
  BOOST_CHECK_THROW(r.at(-1), std::out_of_range);
  BOOST_CHECK_THROW(r.at(3), std::out_of_range);

  r.at(0) = 42;
  r.at(1) = 52;
  r.at(2) = 62;

  BOOST_CHECK_EQUAL(r[0], 42);
  BOOST_CHECK_EQUAL(r[1], 52);
  BOOST_CHECK_EQUAL(r[2], 62);
}

BOOST_AUTO_TEST_CASE(constElementAccess) {
  Int3D a(1, 2,3);
  const Int3DRef r(a);

  BOOST_CHECK_EQUAL(r[0], 1);
  BOOST_CHECK_EQUAL(r[1], 2);
  BOOST_CHECK_EQUAL(r[2], 3);

  BOOST_CHECK_EQUAL(r.at(0), 1);
  BOOST_CHECK_EQUAL(r.at(1), 2);
  BOOST_CHECK_EQUAL(r.at(2), 3);
  BOOST_CHECK_THROW(r.at(-1), std::out_of_range);
  BOOST_CHECK_THROW(r.at(3), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(unaryOps) {
  Int3D a(1, 2, 3);
  Int3D b(42, 52, 62);
  Int3DRef r(a);
  Int3DRef s(b);

  r += s;
  BOOST_CHECK_EQUAL(r[0], 43);
  BOOST_CHECK_EQUAL(r[1], 54);
  BOOST_CHECK_EQUAL(r[2], 65);

  r -= s;
  BOOST_CHECK_EQUAL(r[0], 1);
  BOOST_CHECK_EQUAL(r[1], 2);
  BOOST_CHECK_EQUAL(r[2], 3);

  ConstInt3DRef cs(b);
  r += cs;
  BOOST_CHECK_EQUAL(r[0], 43);
  BOOST_CHECK_EQUAL(r[1], 54);
  BOOST_CHECK_EQUAL(r[2], 65);

  r -= cs;
  BOOST_CHECK_EQUAL(r[0], 1);
  BOOST_CHECK_EQUAL(r[1], 2);
  BOOST_CHECK_EQUAL(r[2], 3);

}

BOOST_AUTO_TEST_CASE(boolOps) {
  Int3D a(1, 2, 3);
  Int3D b(42, 52, 62);
  Int3D c(1, 2, 3);

  Int3DRef r(a);
  Int3DRef s(b);
  Int3DRef t(c);

  BOOST_CHECK_EQUAL(r, r);
  BOOST_CHECK(!(r!=r));

  BOOST_CHECK_NE(r, s);
  BOOST_CHECK(!(r==s));

  BOOST_CHECK_EQUAL(r, t);
  BOOST_CHECK(!(r!=t));

  ConstInt3DRef cr(a);
  ConstInt3DRef cs(b);
  ConstInt3DRef ct(c);
  
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
