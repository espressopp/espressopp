//  (C) Copyright Thomas Brandes, SCAI Fraunhofer

// Boost.Test

// Attention: #include <boost/test/unit_test.hpp> does not include main program 

#include <boost/test/included/unit_test.hpp>

using boost::unit_test_framework::test_suite;

#include "BasicPropertyTest.cpp"

test_suite*
init_unit_test_suite( int argc, char * argv[] ) {

    test_suite* test= BOOST_TEST_SUITE( "Basic Property Test" );

    test->add( new BasicPropertyTestSuite());

    return test;
}

// EOF
