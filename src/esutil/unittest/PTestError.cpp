#define PARALLEL_TEST_MODULE ERROR
#include "include/ut.hpp"

#include "mpi.hpp"
#include "../Error.hpp"

using namespace espresso;
using namespace espresso::esutil;

// Check constructor and no execption at start

BOOST_AUTO_TEST_CASE(constructor) 
{
  Error myError = Error(mpiWorld);

  myError.checkException();
}

static void hangUp()
{
 Error myError = Error(mpiWorld);

 BOOST_TEST_MESSAGE("set exception");

 myError.setException("Hanging exception");
}

// Check for single error

BOOST_AUTO_TEST_CASE(single) 
{
  Error myError = Error(mpiWorld);

  if (mpiWorld->rank() == 0) {

    myError.setException("Set exception");
 
  }

  BOOST_CHECK_THROW(myError.checkException(), 
                    std::runtime_error);
}

// Check set and execption at end 

BOOST_AUTO_TEST_CASE(destructor) 
{
  BOOST_CHECK_THROW(hangUp(), std::runtime_error);
}

