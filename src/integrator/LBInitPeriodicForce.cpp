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
#include "LBInitPeriodicForce.hpp"

namespace espressopp {
  namespace integrator {
    LBInitPeriodicForce::LBInitPeriodicForce(shared_ptr<System> system,
                                             shared_ptr< LatticeBoltzmann > latticeboltzmann)
    : LBInit(system, latticeboltzmann) {
    }

    /* do nothing with initial populations */
    void LBInitPeriodicForce::createDenVel (real _rho0, Real3D u0) {}

    /* set external forces */
    void LBInitPeriodicForce::setForce( Real3D _force ) {
      int _printId = 0;
      int _offset = latticeboltzmann->getHaloSkin();
      Int3D _Ni = latticeboltzmann->getNi();                // system size in lattice node
      Int3D _myNi = latticeboltzmann->getMyNi();          // my local nodes
      Int3D _globIdx = latticeboltzmann->findGlobIdx();   // first lb site global index

      real ampFz = _force[2]; // amplitude of sin force in z-direction
      real perCompX;          // periodicity function with x-dir dependence
      Real3D perF;            // gravity forces in x- and y-dir; in z: amplitude
                              // of force periodic with sin in x, i.e.
                              // perForce_z = ampFz * sin ( x );

      // loop over real lattice nodes //
      real coef = 2. * M_PI / _Ni[0];
      for (int i = _offset; i < _myNi[0] - _offset; i++) {
        perCompX = sin( coef * (_globIdx[0] + i - _offset) );
        perF = Real3D (_force[0], _force[1], ampFz * perCompX);

        for (int j = _offset; j < _myNi[1] - _offset; j++) {
          for (int k = _offset; k < _myNi[2] - _offset; k++) {
            latticeboltzmann->setExtForceLoc( Int3D(i,j,k), perF );
          }
        }
      }

      // take care of the force flag for lb //
      if (_force.sqr() > ROUND_ERROR_PREC) {
        _printId = 1;
        latticeboltzmann->setDoExtForce( true );
        printForce(_force, _printId );
      }
    }

    void LBInitPeriodicForce::addForce(Real3D _force) {
      int _printId = 0;
      int _offset = latticeboltzmann->getHaloSkin();
      Int3D _Ni = latticeboltzmann->getNi();				// system size in lattice node
      Int3D _myNi = latticeboltzmann->getMyNi();		// my local nodes
      Int3D _globIdx = latticeboltzmann->findGlobIdx();// first lb site global index

      bool fOnSite = false, fTot = false;
      real ampFz = _force[2]; // amplitude of sin force in z-direction
      real perCompX;          // periodicity function with x-dir dependence
      Real3D newF, totF;      // gravity forces in x- and y-dir; in z: amplitude
                              // of force periodic with sin in x, i.e.
                              // perForce_z = ampFz * sin ( x );

      // loop over real lattice nodes //
      real coef = 2. * M_PI / _Ni[0];
      for (int i = _offset; i < _myNi[0] - _offset; i++) {
        perCompX = sin( coef * (_globIdx[0] + i - _offset) );
        newF = Real3D (_force[0], _force[1], ampFz * perCompX);

        for (int j = _offset; j < _myNi[1] - _offset; j++) {
          for (int k = _offset; k < _myNi[2] - _offset; k++) {
            totF = latticeboltzmann->getExtForceLoc( Int3D(i,j,k) ) + newF;

            // check if total force is greater than zero
            if (totF.sqr() >  ROUND_ERROR_PREC) {
              fOnSite = true;
              _printId = 2;
              latticeboltzmann->setExtForceLoc( Int3D(i,j,k), totF );
            } else {
              latticeboltzmann->setExtForceLoc( Int3D(i,j,k), Real3D(0.) );
            }

            fTot = fTot | fOnSite; // if at least one site has force -> fTot true
          }
        }
      }

      // take care of the force flag for lb //
      latticeboltzmann->setDoExtForce( fTot );
      printForce(_force, _printId );
    }

    void LBInitPeriodicForce::printForce(Real3D _force, int _printId)
    {
            if (mpiWorld->rank() == 0) {
            using namespace std;

            cout << "External force has been ";
                if (_printId == 0)
              cout << "cancelled. It is now zero.\n" ;
                else if (_printId == 1)
              cout << "set." << endl << "It is a harmonic force with amplitude: "
              << setprecision(5) << _force << endl;
                else
              cout << "added." << endl << "It is a harmonic force with amplitude: "
              << setprecision(5) << _force << endl;

            cout << "-------------------------------------" << endl;

            } else {
                // do nothing
            }
    }

    void LBInitPeriodicForce::registerPython() {
      using namespace espressopp::python;

      class_<LBInitPeriodicForce, bases< LBInit > >
          ("integrator_LBInit_PeriodicForce",	init< shared_ptr< System >,
                                                                                        shared_ptr< LatticeBoltzmann > >())
          .def("setForce", &LBInitPeriodicForce::setForce)
          .def("addForce", &LBInitPeriodicForce::addForce)
      ;
    }
  }
}
