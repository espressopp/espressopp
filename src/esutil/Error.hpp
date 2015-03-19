/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _ESUTIL_ERROR_HPP
#define _ESUTIL_ERROR_HPP

#include <execinfo.h>
#include <stdio.h>
#include "mpi.hpp"

namespace espressopp {
  namespace esutil {

     /** Print Stack*/
    static inline void printStackTrace(std::stringstream &msg, unsigned int max_frames) {

      void *addrlist[max_frames+1];
      int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void *));
      if (addrlen == 0) {
        return;
      }
      msg << "Stack trace:" << std::endl;

      char **symbollist = backtrace_symbols(addrlist, addrlen);
      for(int i = 4; i < addrlen; i++) {
        msg << symbollist[i] << std::endl;
      }
      free(symbollist);
    }

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

        - A processor continues execution if it has set an exception
        - check for exceptions must be invoked explicitly.

    */

    class Error {

    public:

      /** Constructor of the error handler object.

          \param comm is the communicator 

          All processor of the communicator must call
          this constrcutor.
      */

      Error(boost::shared_ptr< boost::mpi::communicator > comm);

      /** Destructor; will also test for pending exceptions. */

      ~Error();

      /** Set an exception. This routine might be called
          individually by one or more processors.

          The error message will only be added if no
          exceptions has been set so far or if add is set
          to true.
      */

      template<typename T>
      void setException(const T msg, bool add = true)
      {
        std::stringstream numMSG;
        numMSG << (noExceptions+1) << "). ";
        if (noExceptions == 0 || add) {
          exceptionMessage += numMSG.str();
          exceptionMessage += msg;
          exceptionMessage += "\n";
        } 
        noExceptions++;
      }

      /** Checks whether any processor has set an exception and
          if yes will throw a C++ exception.

          The string of the exception will be the concatenation
          of all messages set so far.
      */

      void checkException();

    private:    

      boost::shared_ptr< boost::mpi::communicator > comm;

      std::string exceptionMessage;

      int noExceptions;  //!< counts exceptions on this proc
    };
    
    
  }
}

#endif
