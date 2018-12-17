/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
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

// ESPP_CLASS
#ifndef _INTEGRATOR_MDINTEGRATOR_HPP
#define _INTEGRATOR_MDINTEGRATOR_HPP
#include "python.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"
#include "Extension.hpp"
#include <boost/signals2.hpp>
#include "types.hpp"
#include "esutil/Error.hpp"


namespace espressopp {
  namespace integrator {

    /** Abstract base class for Molecular Dynamics Integrator.

        Note: Class accesses system object for storage, bc, communicator,
              interaction list.
    */

    class Extension; //fwd declaration

    struct ExtensionList : public std::vector<shared_ptr<Extension> > {
          typedef esutil::ESPPIterator<std::vector<Extension> > Iterator;
    };

    class MDIntegrator : public SystemAccess {
      public:
        /** Constructor for an integrator.
            \param system is reference to the system.
            Note: This class will keep a weak reference to the system.
        */
        MDIntegrator(shared_ptr<System> system);

        /** Destructor. */
        virtual ~MDIntegrator();

        /** Setter routine for the timestep. */
        void setTimeStep(real dt);

        /** Getter routine for the timestep. */
        real getTimeStep() { return dt; }

        /** Getter routine for integration step */
        void setStep(long long step_) { step = step_; }

        /** Getter routine for integration step */
        long long getStep() { return step; }

        /** This method runs the integration for a certain number of steps. */
        virtual void run(int nsteps) = 0;

        void addExtension(shared_ptr<integrator::Extension> extension);

        shared_ptr<integrator::Extension> getExtension(int k);

        int getNumberOfExtensions();

        // signals to extend the integrator
        boost::signals2::signal<void ()> runInit; // initialization of run()
        boost::signals2::signal<void ()> recalc1; // inside recalc, before updateForces()
        boost::signals2::signal<void ()> aftCalcSlow; // after calculation of slow forces updateForces(true) in VerlocityVerletRESPA
        boost::signals2::signal<void ()> recalc2; // inside recalc, after  updateForces()
        boost::signals2::signal<void ()> befIntP; // before integrate1()
        boost::signals2::signal<void (real&)> inIntP; // inside end of integrate1()
        boost::signals2::signal<void ()> aftIntP; // after  integrate1()
        boost::signals2::signal<void ()> aftInitF; // after initForces()
        boost::signals2::signal<void ()> aftCalcF; // after calcForces()
        boost::signals2::signal<void ()> befIntV; // before integrate2()
        boost::signals2::signal<void ()> aftIntV; // after  integrate2()
        boost::signals2::signal<void ()> aftIntSlow; // after integrateSlow() in VerlocityVerletRESPA


        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        bool timeFlag;

        ExtensionList exList;

        /** Integration step */
        long long step;

        /** Timestep used for integration */
        real dt;

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };


  }
}

#endif
