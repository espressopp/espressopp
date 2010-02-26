#include "Error.hpp"

#include <stdexcept>
#include <cstdlib>
#include <sstream>

using namespace espresso;
using namespace esutil;
using namespace boost;

/********************************************************************/

Error::Error(mpi::communicator _comm)
{
  comm = _comm;
  noExceptions = 0;
}

/********************************************************************/

Error::~Error()
{
  // final test for any hanging exception

  checkException();
}

/********************************************************************/

void Error::setException(const char* msg, bool add)
{
  // do not add message if not first and add == false

  if (noExceptions == 0 || add) {
    exceptionMessage += msg;
  } 

  noExceptions++;
}

/********************************************************************/

void Error::checkException()
{
  // count exceptions on all processors

  int totalExceptions = 0;

  mpi::all_reduce(comm, noExceptions, totalExceptions, std::plus<int>());

  printf("On proc %d: noExceptions = %d, total = %d\n",
          comm.rank(), noExceptions, totalExceptions);

  if (totalExceptions > 0) {

    std::ostringstream msg;

    msg << totalExceptions << " exceptions occurred";

    if (exceptionMessage.length() > 0) {
      msg << ": " << exceptionMessage;
      exceptionMessage.clear();
    }

    // reset exception counter for next time

    noExceptions = 0;

    throw std::runtime_error(msg.str());
  }
}

