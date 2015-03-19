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
#include "NPart.hpp"
#include "storage/DomainDecomposition.hpp"

using namespace espressopp;

namespace espressopp {
  namespace analysis {
    real NPart::compute_real() const {

      int myN, systemN;
      System& system = getSystemRef();
      myN = system.storage->getNRealParticles();
      myN += system.storage->getNAdressParticles();
      boost::mpi::reduce(*getSystem()->comm, myN, systemN, std::plus<int>(), 0);
      
      return 1.0*systemN;
    }

    void NPart::registerPython() {
      using namespace espressopp::python;
      class_<NPart, bases< Observable > >
        ("analysis_NPart", init< shared_ptr< System > >())
      ;
    }
  }
}
