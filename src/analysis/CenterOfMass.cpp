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
#include <cmath>
#include "CenterOfMass.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espressopp;
using namespace iterator;

namespace espressopp {
  namespace analysis {
    
    // ! it's folded COM
    Real3D CenterOfMass::computeVector() const {

      System& system = getSystemRef();
  
      // compute the center-of-mass of the system (r_CM = \sum m_i r_i / \sum m_i)
      real xcom = 0.0;
      real ycom = 0.0;
      real zcom = 0.0;
      real mass = 0.0;
      real xcom_sum = 0.0;
      real ycom_sum = 0.0;
      real zcom_sum = 0.0;
      real mass_sum = 0.0;

      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        const Particle& p = *cit;
        xcom += p.mass() * p.position()[0];
        ycom += p.mass() * p.position()[1];
        zcom += p.mass() * p.position()[2];
        mass += p.mass();
      }
      
      // it was reduce, but we need it for all cpus
      boost::mpi::all_reduce(*mpiWorld, xcom, xcom_sum, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, ycom, ycom_sum, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, zcom, zcom_sum, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, mass, mass_sum, std::plus<real>());

      //Real3D force(0.0, 0.0, 0.0);
      return Real3D(xcom_sum / mass_sum, ycom_sum / mass_sum, zcom_sum / mass_sum);
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real CenterOfMass::compute() const {
      return -1.0;
    }

    void CenterOfMass::registerPython() {
      using namespace espressopp::python;
      class_<CenterOfMass, bases< Observable > >
        ("analysis_CenterOfMass", init< shared_ptr< System > >())
        .def("compute", &CenterOfMass::computeVector)
      ;
    }
  }
}
