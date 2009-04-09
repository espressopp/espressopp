#include <acconfig.hpp>
#define BOOST_TEST_MODULE PropertyVector

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/foreach.hpp>

#include "../TupleVector.hpp"

using namespace espresso::esutil;

struct Fixture {
    TupleVector mv;
    const TupleVector &constMv;

    TupleVector::PropertyId intProp, floatProp;

    Fixture(): constMv(mv) {
	intProp = mv.addProperty<int>();
	floatProp = mv.addProperty<float>(3);
        mv.setGranularity(8);
        mv.setShrinkThreshold(32);
	mv.resize(42);
    }

    ~Fixture() {
    }
};

//____________________________________________________________________________//


BOOST_FIXTURE_TEST_CASE(properties_resize_test, Fixture)
{
    // check that empty vector is empty
    TupleVector mve;
    mve.setGranularity(8);
    mve.setShrinkThreshold(32);
    BOOST_CHECK(mve.size() == 0);
    BOOST_CHECK(mve.getNumProperties() == 0);

    // but not, if we add a Property
    mve.addProperty<int>(3);
    BOOST_REQUIRE(mve.getNumProperties() == 1);

    // add a second property
    TupleVector::PropertyId prop = mve.addProperty<float>();
    BOOST_REQUIRE(mve.getNumProperties() == 2);

    // now, delete the property again
    mve.eraseProperty(prop);
    BOOST_REQUIRE(mve.getNumProperties() == 1);

    // check resizing copy
    TupleVector mvc(mve, 42);
    BOOST_CHECK_EQUAL(mvc.size(), size_t(42));
    BOOST_CHECK_EQUAL(mvc.capacity(), size_t(48));
    BOOST_CHECK_EQUAL(mvc.getNumProperties(), size_t(1));
}

BOOST_FIXTURE_TEST_CASE(particles_resize_test, Fixture)
{
    // this test needs the initial size to be 42
    BOOST_REQUIRE(mv.size() == 42);

    // resizing above threshold, should shrink
    mv.clear();
    BOOST_CHECK(mv.size() == 0);
    BOOST_CHECK(mv.capacity() == 0);

    mv.resize(12);
    BOOST_CHECK(mv.size() == 12);
    BOOST_CHECK(mv.capacity() == 16);

    // resizing below the threshold should not shrink
    mv.clear();
    BOOST_CHECK(mv.size() == 0);
    BOOST_CHECK(mv.capacity() == 16);
}

BOOST_FIXTURE_TEST_CASE(pointers_test, Fixture)
{
    TupleVector::PropertyReference<int> pRef = mv.getProperty<int>(intProp);

    ////// element references
    TupleVector::reference ref = mv[17];
    TupleVector::const_reference constRef = constMv[17];
    TupleVector::pointer ptr = &ref;
    TupleVector::const_pointer constPtr = &constRef;
    // check comparisons
    BOOST_CHECK(ptr == constPtr);
    BOOST_CHECK(!(ptr != constPtr));

    // check changing the pointer
    constPtr = &mv[3];
    BOOST_CHECK(ptr != constPtr);

    // this does not compile, tries const -> non-const conversion
    //TupleVector::pointer ptr2 = &constRef;

    // check that addresses seem to match
    pRef[ref] = 42;
    pRef[*ptr] = 45;    
    BOOST_CHECK(pRef[ref] == 45);
}

BOOST_FIXTURE_TEST_CASE(dereference_scalar_test, Fixture)
{
    ////// property references
    const TupleVector &constMv = mv;
    TupleVector::PropertyReference<int> pRef = mv.getProperty<int>(intProp);
    TupleVector::ConstPropertyReference<int> constPRef = constMv.getProperty<int>(intProp);

    // convert non-const -> const
    { TupleVector::ConstPropertyReference<int> constPRef2 = pRef; }

    // this does not compile, tries const -> non-const conversion
    //TupleVector::PropertyReference<int> pRef2 = constMv.getProperty<int>(intProp);    

    ////// element references
    TupleVector::reference ref = mv[0];
    TupleVector::const_reference constRef = constMv[0];

    // this does not compile, reassigning a non-trivial reference
    // ref = ref;

    // convert non-const -> const
    { TupleVector::const_reference constRef = mv[0]; }

    // this does not compile, tries const -> non-const conversion
    //TupleVector::reference ref2 = constMv[0];

    BOOST_CHECK_THROW(TupleVector::reference ref2 = mv.at(42),
		      std::out_of_range);

    BOOST_CHECK_THROW(TupleVector::const_reference ref2 = constMv.at(43),
		      std::out_of_range);

    pRef[ref] = 176;
    BOOST_CHECK_EQUAL(pRef[ref], 176);
    pRef[ref] = 42;
    BOOST_CHECK_EQUAL(constPRef[constRef], 42);

    // this does not compile, overriding const in various ways
    //pRef[constRef] = 42;
    //constPRef[ref] = 42;
    //constPRef[constRef] = 42;
    //TupleVector::ReferenceIndexAccess::getIndex(constRef);
}

BOOST_FIXTURE_TEST_CASE(dereference_array_test, Fixture)
{
    ////// property references
    const TupleVector &constMv = mv;
    TupleVector::ArrayPropertyReference<float>           pRef = mv.getArrayProperty<float>(floatProp);
    TupleVector::ConstArrayPropertyReference<float> constPRef = constMv.getArrayProperty<float>(floatProp);

    // convert non-const -> const
    { TupleVector::ConstArrayPropertyReference<float> constPRef2 = pRef; }

    // this does not compile, tries const -> non-const conversion
    //TupleVector::ArrayPropertyReference<float, 3> pRef2 = constMv.getArrayProperty<float, 3>(floatProp);

    ////// element references
    TupleVector::reference ref = mv[0];
    TupleVector::const_reference constRef = constMv[0];

    // this does not compile, tries const -> non-const conversion
    //TupleVector::reference ref2 = constMv[0];

    BOOST_CHECK_THROW(TupleVector::reference ref2 = mv.at(42),
		      std::out_of_range);

    BOOST_CHECK_THROW(TupleVector::const_reference ref2 = constMv.at(43),
		      std::out_of_range);

    // to check that array elements do not overlap, create two adjacent elements
    TupleVector::reference ref2 = mv[1];

    pRef[ref][0] = 0.01;
    pRef[ref][1] = 0.2;
    pRef[ref][2] = 3.0;

    pRef[ref2][0] = 0.04;
    pRef[ref2][1] = 0.5;
    pRef[ref2][2] = 6.0;

    BOOST_CHECK_CLOSE(constPRef[constRef][0], 0.01f, 1e-10f);
    BOOST_CHECK_CLOSE(          pRef[ref][1], 0.20f, 1e-10f);
    BOOST_CHECK_CLOSE(     pRef[constRef][2], 3.00f, 1e-10f);

    BOOST_CHECK_CLOSE(constPRef[ref2][0], 0.04f, 1e-10f);
    BOOST_CHECK_CLOSE(constPRef[ref2][1], 0.50f, 1e-10f);
    BOOST_CHECK_CLOSE(constPRef[ref2][2], 6.00f, 1e-10f);

    // this does not compile, overriding const in various ways
    //pRef[constRef][0] = 42;
    //constPRef[ref][1] = 42;
    //constPRef[constRef][2] = 42;
}

BOOST_FIXTURE_TEST_CASE(iterator_test, Fixture)
{
  {
    TupleVector::PropertyReference<int> intRef = mv.getProperty<int>(intProp);

    // fill int property
    size_t i = 0;
    BOOST_FOREACH(TupleVector::reference ref, mv) {
	intRef[ref] = i++;
    }
    BOOST_CHECK_EQUAL(i, mv.size());

    // and check the data is there, with constant array this time
    i = 0;
    BOOST_FOREACH(TupleVector::const_reference ref, constMv) {
	BOOST_CHECK_EQUAL(intRef[ref], int(i++));
    }
    // this code above should not compile without const_
  }

  // insert two elements by insert()
  TupleVector::iterator it = mv.begin() + 10;
  // iterator it actually never changes, because it points to the
  // inserted element
  it = mv.insert(it);
  BOOST_CHECK_EQUAL(size_t(43), mv.size());
  it = mv.insert(it, mv[6]);
  BOOST_CHECK_EQUAL(size_t(44), mv.size());
  it = mv.insert(it, mv[3]);
  BOOST_CHECK_EQUAL(size_t(45), mv.size());
  mv.insert(it,3);
  BOOST_CHECK_EQUAL(size_t(48), mv.size());

  /* check the data is there, with constant array this time.
     After above magic, we have:
     mv[0:9]=0:9
     mv[10:12] = x
     mv[13] = mv[3] = 3
     mv[14] = mv[6] = 6
     mv[15] = x
     mv[16:] = 10:
  */
  {
    TupleVector::PropertyReference<int> intRef = mv.getProperty<int>(intProp);
    size_t i = 0;
    BOOST_FOREACH(TupleVector::const_reference ref,
		  static_cast<const TupleVector &>(mv)) {
	if (i < 10) {
	    BOOST_CHECK_EQUAL(intRef[ref], int(i));
	} else if (i == 13) {
	    BOOST_CHECK_EQUAL(intRef[ref], 3);
	} else if (i == 14) {
	    BOOST_CHECK_EQUAL(intRef[ref], 6);
	} else if (i >= 16) {
	    BOOST_CHECK_EQUAL(intRef[ref], int(i - 6));
	}
        i++;
    }
  }

  // and now, delete some elements
  it = mv.begin() + 10;
  it = mv.erase(it, it + 3);
  BOOST_CHECK_EQUAL(size_t(45), mv.size());
  it = mv.erase(it + 2);
  BOOST_CHECK_EQUAL(size_t(44), mv.size());

  /* check the data is there, with constant array this time.
     Now, we have:
     mv[0:9]=0:9
     mv[10] = mv[3] = 3
     mv[11] = mv[6] = 6
     mv[12:] = 10:
  */
  {
    TupleVector::PropertyReference<int> intRef = mv.getProperty<int>(intProp);

    size_t i = 0;
    BOOST_FOREACH(TupleVector::const_reference ref,
		  static_cast<const TupleVector &>(mv)) {
	if (i < 10) {
	    BOOST_CHECK_EQUAL(intRef[ref], int(i));
	} else if (i == 10) {
	    BOOST_CHECK_EQUAL(intRef[ref], 3);
	} else if (i == 11) {
	    BOOST_CHECK_EQUAL(intRef[ref], 6);
	} else if (i >= 12) {
	    BOOST_CHECK_EQUAL(intRef[ref], int(i - 2));
	}
        i++;
    }
  }

  // and now, move some elements
  mv.copy(mv[10], mv[9]);
  BOOST_CHECK_EQUAL(size_t(44), mv.size());
  mv.copy(mv.begin() + 3, mv.begin() + 7, mv.begin() + 11);
  BOOST_CHECK_EQUAL(size_t(44), mv.size());

  /* check the data is there, with constant array this time.
     Now, we have:
     mv[0:9]=0:9
     mv[10] = mv[9] = 9
     mv[11:14] = mv[3:6] = 3:6
     mv[15:] = 13:
  */
  {
    TupleVector::PropertyReference<int> intRef = mv.getProperty<int>(intProp);

    size_t i = 0;
    BOOST_FOREACH(TupleVector::const_reference ref,
		  static_cast<const TupleVector &>(mv)) {
      if (i < 10) {
	BOOST_CHECK_EQUAL(intRef[ref], int(i));
      } else if (i == 10) {
	BOOST_CHECK_EQUAL(intRef[ref], 9);
      } else if (i < 15) {
	BOOST_CHECK_EQUAL(intRef[ref], int(i - 8));
      } else {
	BOOST_CHECK_EQUAL(intRef[ref], int(i - 2));
      }
      i++;
    }
  }
}

BOOST_FIXTURE_TEST_CASE(allocation_test, Fixture)
{
  mv.clear();

  // stress-test by continous growing
  // checks realloc code
  for (size_t i = 0; i < 10000; ++i)
    mv.insert(mv.end());

  BOOST_CHECK_EQUAL(size_t(10000), mv.size());

  // add property to new field
  // checks malloc code
  TupleVector::PropertyId newProp = mv.addProperty<float>(3);

  // try wether writing to present and new fields raises memory problems
  {
    TupleVector::PropertyReference<int> intRef = mv.getProperty<int>(intProp);
    TupleVector::ArrayPropertyReference<float> fltRef
      = mv.getArrayProperty<float>(floatProp);
    TupleVector::ArrayPropertyReference<float> fltNRef
      = mv.getArrayProperty<float>(newProp);

    
    BOOST_FOREACH(TupleVector::reference ref, mv) {
      intRef[ref] = 42;
      fltRef[ref][0] =   4.2;
      fltRef[ref][1] =  42.0;
      fltRef[ref][2] = 420.0;
      fltNRef[ref][0] =   1.2;
      fltNRef[ref][1] =  12.0;
      fltNRef[ref][2] = 120.0;
    }

    BOOST_FOREACH(TupleVector::reference ref, mv) {
      BOOST_CHECK_EQUAL(intRef[ref], 42);
      BOOST_CHECK_CLOSE(fltRef[ref][0],   4.2f, 1e-4f);
      BOOST_CHECK_CLOSE(fltRef[ref][1],  42.0f, 1e-4f);
      BOOST_CHECK_CLOSE(fltRef[ref][2], 420.0f, 1e-4f);
      BOOST_CHECK_CLOSE(fltNRef[ref][0],   1.2f, 1e-4f);
      BOOST_CHECK_CLOSE(fltNRef[ref][1],  12.0f, 1e-4f);
      BOOST_CHECK_CLOSE(fltNRef[ref][2], 120.0f, 1e-4f);
    }
  }

  mv.eraseProperty(newProp);

  for (size_t i = 0; i < 10000; ++i)
    mv.erase(mv.begin());

  BOOST_CHECK_EQUAL(size_t(0), mv.size());
}
