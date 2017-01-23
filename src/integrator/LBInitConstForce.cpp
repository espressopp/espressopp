/*
  Copyright (C) 2012-2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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
#include <iomanip>
#include "LBInitConstForce.hpp"

namespace espressopp {
  namespace integrator {
//    LOG4ESPP_LOGGER(LBInitConstForce::theLogger, "LBInitConstForce");
    LBInitConstForce::LBInitConstForce(shared_ptr<System> system,
                                       shared_ptr< LatticeBoltzmann > latticeboltzmann)
    : LBInit(system, latticeboltzmann) {
    }

    /* do nothing with initial populations, as we deal here only with forces */
    void LBInitConstForce::createDenVel (real _rho0, Real3D u0) {}

		/* SET BODY FORCE */
    void LBInitConstForce::setForce(Real3D _force)
    {
      int _id = 0;		// _id shows whether force had been cancelled, set or added
			int _offset = latticeboltzmann->getHaloSkin();
      Int3D _Ni = latticeboltzmann->getMyNi();

      /* external forces loop */
      for (int i = _offset; i < _Ni.getItem(0)-_offset; i++) {
        for (int j = _offset; j < _Ni.getItem(1)-_offset; j++) {
          for (int k = _offset; k < _Ni.getItem(2)-_offset; k++) {
            // set local forces and general flag
            if (_force != Real3D(0.,0.,0.)) {
              latticeboltzmann->setDoExtForce(true);
              latticeboltzmann->setExtForceLoc(Int3D(i,j,k),_force);
              _id = 1;
            } else {
              latticeboltzmann->setDoExtForce(false);
              latticeboltzmann->setExtForceLoc(Int3D(i,j,k),Real3D(0.,0.,0.));
            }
          }
        }
      }
      printForce(_force, _id);
    }

		/* ADD BODY FORCE TO EXISTING FORCES (IF ANY) */
    void LBInitConstForce::addForce(Real3D _force)
    {
      int _id = 0;
			int _offset = latticeboltzmann->getHaloSkin();
			Int3D _Ni = latticeboltzmann->getMyNi();

      Real3D existingforce;

      /* external forces loop */
			for (int i = _offset; i < _Ni.getItem(0)-_offset; i++) {
				for (int j = _offset; j < _Ni.getItem(1)-_offset; j++) {
					for (int k = _offset; k < _Ni.getItem(2)-_offset; k++) {
            existingforce = latticeboltzmann->getExtForceLoc(Int3D(i,j,k));
            // set local forces and general flag
            if (existingforce + _force != Real3D(0.,0.,0.)) {
              latticeboltzmann->setDoExtForce(true);
              latticeboltzmann->setExtForceLoc(Int3D(i,j,k),existingforce + _force);
              _id = 2;
            } else {
              latticeboltzmann->setDoExtForce(false);
              latticeboltzmann->setExtForceLoc(Int3D(i,j,k),Real3D(0.,0.,0.));
            }
          }
        }
      }
      printForce(_force, _id);
    }

    void LBInitConstForce::printForce(Real3D _force, int _id)
    {
			if (mpiWorld->rank()==0) {
				// print constant force
				using std::setprecision;
				using std::fixed;
				using std::setw;
				
				std::cout << setprecision(5);
				std::cout << "-------------------------------------\n";
				
				if (_id == 0) {
					std::cout << "External force has been cancelled. It is now zero.\n" ;
				} else if (_id == 1) {
					std::cout << "External force has been set. It is a constant force:\n" ;
				} else if (_id == 2) {
					std::cout << "External force has been added. It is a constant force:\n" ;
				} else {
				}
				
				if (_id != 0) {
					std::cout << " extForce.x is " << _force.getItem(0) << "\n";
					std::cout << " extForce.y is " << _force.getItem(1) << "\n";
					std::cout << " extForce.z is " << _force.getItem(2) << "\n";
					std::cout << "-------------------------------------\n";
				}
			} else {
				// do nothing
			}
    }

    void LBInitConstForce::registerPython() {
      using namespace espressopp::python;

      class_<LBInitConstForce, bases< LBInit > >
          ("integrator_LBInit_ConstForce", init<	shared_ptr< System >,
																							shared_ptr< LatticeBoltzmann > >())
          .def("setForce", &LBInitConstForce::setForce)
          .def("addForce", &LBInitConstForce::addForce)
      ;
    }
  }
}

