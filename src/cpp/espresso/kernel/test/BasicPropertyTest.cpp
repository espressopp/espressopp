//  (C) Copyright Thomas Brandes, SCAI Fraunhofer

// Boost.Test

// Attention: #include <boost/test/unit_test.hpp> does not include main program 

#include <boost/test/included/unit_test.hpp>

using boost::unit_test_framework::test_suite;
using boost::unit_test_framework::test_case;

// the class to be tested

#include "../BasicProperty.hpp"

// needed interface to create mock object of ParticleContainer

#include "../ParticleContainer.hpp"

#include <iostream>

using namespace std;

class BasicPropertyTest {

  private: 

    BasicProperty *p;

  public:

  // constructor

  BasicPropertyTest() {

    p = new ArrayProperty<double,3> ("velocity");

  }

  /** test attributes  */

  void test1() {

    /* critical tests for attributes */

    BOOST_CHECK(p->getSize() == sizeof(double) * 3); // critical test 
    BOOST_CHECK(p->getName() == "velocity"); // critical test 

  }

  /** test registration  */

  void test2() {

    /* test for registration of properties in a particle container */

    ParticleContainer *c = new ParticleContainer (0, 0);

    p->registerContainer (c, 2);
    BOOST_CHECK(p->isRegistered(c));
    BOOST_CHECK(p->getPosition() == 2);
    p->updatePosition(3);
    BOOST_CHECK(p->getPosition() == 3);

    // this should result in an ASSERTION error 
    BOOST_CHECK_THROW(p->registerContainer (c, 2), std::bad_exception)
    p->deregister ();
    BOOST_CHECK(p->isRegistered(c) == false);
  }

  /** test data access */

  void test3() {

    ParticleContainer *c = new ParticleContainer (0, 1);
    p->registerContainer (c, 0);
    ParticleRef particle = c->addParticle(0, 0);
    double vals [3] = { 1.0, 1.1, 1.2 };
    double (*hVals) [3] =  (double (*) [3]) p->getPtrData (particle);
    double (&pVals) [3] = hVals[0];

    p->storeData (particle, vals);
    BOOST_CHECK_EQUAL (pVals[0], 1.0);
    pVals [1] = 1.3;
    BOOST_CHECK_NE (pVals[1], vals[1]);
    p->fetchData (particle, vals); 
    BOOST_CHECK_EQUAL (vals[1], pVals[1]);
    BOOST_CHECK_EQUAL (vals[2], pVals[2]);
  }

  ~BasicPropertyTest() {

    delete(p);

  }

};

class BasicPropertyTestSuite : public test_suite
{
   public:

   BasicPropertyTestSuite() : test_suite("BasicProperty test suite") {

      // create an instance of the test cases class
      boost::shared_ptr<BasicPropertyTest> instance(new BasicPropertyTest());

      // create the test cases
      test_case* myTestCase1 = BOOST_CLASS_TEST_CASE(
                              &BasicPropertyTest::test1, instance );

      // add the test cases to the test suite

      add(myTestCase1);
 
      add (BOOST_CLASS_TEST_CASE(&BasicPropertyTest::test2, instance ));
      add (BOOST_CLASS_TEST_CASE(&BasicPropertyTest::test3, instance ));

   }

};

test_suite*
init_unit_test_suite( int argc, char * argv[] ) {

    test_suite* test= BOOST_TEST_SUITE( "Basic Property Test" );

    test->add( new BasicPropertyTestSuite());

    return test;
}

// EOF
