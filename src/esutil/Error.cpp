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

#include "Error.hpp"

#include <stdexcept>
#include <cstdlib>
#include <sstream>

using namespace boost;

namespace espressopp {
  namespace esutil {

    /********************************************************************/
    Error::Error(shared_ptr< mpi::communicator > _comm)
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

    void Error::checkException()
    {
      // count exceptions on all processors

      int totalExceptions = 0;

      mpi::all_reduce(*comm, noExceptions, totalExceptions, std::plus<int>());

      //printf("On proc %d: noExceptions = %d, total = %d\n",
	  //   comm->rank(), noExceptions, totalExceptions);

      if (totalExceptions > 0) {

        std::ostringstream msg;

        msg << totalExceptions << " exceptions occurred";

        if (exceptionMessage.length() > 0) {
          msg << ":\n cpu "<< comm->rank()<< ":  Exception message(s):\n" << exceptionMessage;
          msg <<"\n";
          msg << "On proc "<< comm->rank()<< ": exceptions = "<<noExceptions<<
                  ", total = "<< totalExceptions <<"\n";
 
          exceptionMessage.clear();
        }

        // reset exception counter for next time

        noExceptions = 0;

        throw std::runtime_error(msg.str());
      }
    }

  }
}
