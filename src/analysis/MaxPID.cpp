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

#include "python.hpp"
#include "MaxPID.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espressopp;
using namespace iterator;

namespace espressopp {
  namespace analysis {
    real MaxPID::compute() const {

      long myMaxPID, systemMaxPID;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      myMaxPID = 0;
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        long pid = cit->id();
        if (pid > myMaxPID) {
          myMaxPID = pid;
        }
      }

      // it was reduce
      boost::mpi::all_reduce(*getSystem()->comm, myMaxPID, systemMaxPID, mpi::maximum<long>());

      return 1.0*systemMaxPID;

    }

    void MaxPID::registerPython() {
      using namespace espressopp::python;
      class_<MaxPID, bases< Observable > >
        ("analysis_MaxPID", init< shared_ptr< System > >())
      ;
    }
  }
}
