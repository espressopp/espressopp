#ifndef _ERROR_HPP
#define _ERROR_HPP

#include "mpi.hpp"

namespace espresso {
  namespace esutil {


  /** This class supports exception handling for MPI processors.

      The problem of MPI is that it is not possible to throw
      an exception on individual processors because it might
      seriously damage running communication patterns, e.g.
      they will not participate anymore in a global communication
      and the other processors will hang.

      So it is only safe to throw an exception if this is done on
      all processors.

      This class allows processors to throw individual exceptions
      and a global exception will be thrown later.

      Be careful about the drawbacks:

      - A processor continues exectution if it has set an exception
      - check for exceptions must be invoked explicitly.

  */

  class Error {

  public:

    /** Constructor of the error handler object.

        \param comm is the communicator 

        All processor of the communicator must call
        this constrcutor.
    */

    Error(boost::mpi::communicator comm);

    /** Destructor; will also test for pending exceptions. */

    ~Error();

    /** Set an exception. This routine might be called
        individually by one or more processors.

        The error message will only be added if no
        exceptions has been set so far or if add is set
        to true.
    */

    void setException(const char* msg, bool add = true);

    /** Checks whether any processor has set an exception and
        if yes will throw a C++ exception.

        The string of the exception will be the concatenation
        of all messages set so far.
    */

    void checkException();

  private:    

    boost::mpi::communicator comm;

    std::string exceptionMessage;

    int noExceptions;  //!< counts exceptions on this proc
  }; 
}
}

#endif
