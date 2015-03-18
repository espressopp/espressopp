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
#include "VelocityVerletOnRadius.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(VelocityVerletOnRadius::theLogger, "VelocityVerletOnRadius");

    VelocityVerletOnRadius::VelocityVerletOnRadius(shared_ptr<System> system, real _radialDampingMass)
      : Extension(system), radialDampingMass(_radialDampingMass)
    {
      LOG4ESPP_INFO(theLogger, "VelocityVerletOnRadius constructed");
    }

    void VelocityVerletOnRadius::disconnect(){
      _aftIntP.disconnect();
      _aftIntV.disconnect();
      _aftInitF.disconnect();
    }

    void VelocityVerletOnRadius::connect(){
      // connection to initialisation
      _aftIntP  = integrator->aftIntP.connect( boost::bind(&VelocityVerletOnRadius::integrate1, this));
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&VelocityVerletOnRadius::integrate2, this));
      _aftInitF  = integrator->aftInitF.connect( boost::bind(&VelocityVerletOnRadius::initForces, this));
    }

    void VelocityVerletOnRadius::integrate1() {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      // loop over all particles of the local cells
      real dt = integrator->getTimeStep();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", forceRadius = " << cit->fradius());
      }
     for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        //LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() <<  ", radius = " << cit->radius());
        real dtfm = 0.5 * dt / radialDampingMass;
        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
        cit->vradius() += dtfm * cit->fradius();
        // Propagate radius (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
        real deltaP = cit->vradius();
        deltaP *= dt;
        cit->radius() += deltaP;
      }
    }

    void VelocityVerletOnRadius::integrate2() {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      // loop over all particles of the local cells
      real dt = integrator->getTimeStep();
      real half_dt = 0.5 * dt;
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real dtfm = half_dt / radialDampingMass;
        /* Propagate radius velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        cit->vradius() += dtfm * cit->fradius();
      }
    }

    void VelocityVerletOnRadius::initForces() {
      // fradius is initialized for real + ghost particles
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      LOG4ESPP_INFO(theLogger, "init fradius for real particles");
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        cit->fradius() = 0.0;
      }
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletOnRadius::registerPython() {

      using namespace espressopp::python;

      class_<VelocityVerletOnRadius, shared_ptr<VelocityVerletOnRadius>, bases<Extension> >
        ("integrator_VelocityVerletOnRadius", init< shared_ptr< System >, real >())
        .def("connect", &VelocityVerletOnRadius::connect)
        .def("disconnect", &VelocityVerletOnRadius::disconnect)
        .add_property("radialDampingMass" ,&VelocityVerletOnRadius::getRadialDampingMass, &VelocityVerletOnRadius::setRadialDampingMass )
        ;
    }

  }
}
