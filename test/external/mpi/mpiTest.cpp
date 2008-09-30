//  (C) Copyright Thomas Brandes, SCAI Fraunhofer

// MpiTest

// Attention: #include <boost/test/unit_test.hpp> does not include main program

#include <boost/test/included/unit_test.hpp>
#include "mpi.h"

using boost::unit_test_framework::test_suite;

// the class to be tested

void test_init() {

   // make sure that we run on more than 1 processor

   int mpiSize = MPI::COMM_WORLD.Get_size();
   BOOST_CHECK(mpiSize > 1);

}

#define BROADCAST_VAL 13

void test_bcast() {

   int val = 0;
   int rank = MPI::COMM_WORLD.Get_rank();

   if (rank == 0) {
      val = BROADCAST_VAL;
   }

   MPI::COMM_WORLD.Bcast(&val, 1, MPI_INT, 0);

   // verify that now all processors have the value

   BOOST_CHECK_EQUAL (val, BROADCAST_VAL);
}


void test_sum() {

   int val  = MPI::COMM_WORLD.Get_rank() + 1;

   int sum  = 1711;  // must be overwritten

   // val is send buffer, sum is receive buffer 

   MPI::COMM_WORLD.Allreduce(&val, &sum, 1, MPI_INT, MPI_SUM);

   int size = MPI::COMM_WORLD.Get_size();

   size = size * (size + 1) / 2;

   BOOST_CHECK_EQUAL(sum, size);
}

class MPITestSuite : public test_suite
{
   public:

   MPITestSuite() : test_suite("Simple MPI test suite") {

      add (BOOST_TEST_CASE(&test_init));
      add (BOOST_TEST_CASE(&test_bcast));
      add (BOOST_TEST_CASE(&test_sum));

   }

};

test_suite*
init_unit_test_suite( int argc, char * argv[] ) {

    MPI::Init(argc, argv);

    test_suite* test= BOOST_TEST_SUITE( "MPI Test" );

    test->add( new MPITestSuite());

    return test;
}

