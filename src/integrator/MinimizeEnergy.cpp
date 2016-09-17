/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)

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
#include "MinimizeEnergy.hpp"

namespace espressopp {
    namespace integrator {
	
	using namespace interaction;
	using namespace iterator;
	using namespace esutil;
	
	LOG4ESPP_LOGGER(MinimizeEnergy::theLogger, "MinimizeEnergy");
	
	MinimizeEnergy::MinimizeEnergy(shared_ptr<System> system,
				       real gamma,
				       real ftol,
				       real max_displacement,
				       bool variable_step_flag)
	    
	    : SystemAccess(system), gamma_(gamma), ftol_sqr_(ftol),
	      max_displacement_(max_displacement), variable_step_flag_(variable_step_flag) {
	    LOG4ESPP_INFO(theLogger, "construct MinimizeEnergy");
	    resort_flag_ = true;
	    dp_MAX = 0.;
	    nstep_ = 0;
	}
	
	MinimizeEnergy::~MinimizeEnergy() {
	    LOG4ESPP_INFO(theLogger, "free MinimizeEnergy");
	}
	
	bool MinimizeEnergy::run(int max_steps, bool verbose) {
	    bool retval = false;
	    System &system = getSystemRef();
	    storage::Storage &storage = *system.storage;
	    real skin_half = 0.5 * system.getSkin();
	    dp_sqr_max_ = 0.0;
	    f_max_sqr_ = std::numeric_limits<real>::max();
	    
	    // Before start make sure that particles are on the right processor
	    if (resort_flag_) {
		LOG4ESPP_DEBUG(theLogger, "storage.decompose")
		    storage.decompose();
		resort_flag_ = false;
	    }
	    
	    updateForces();
	    
	    LOG4ESPP_INFO(theLogger,
			  "starting energy minimalization loop (iters=" << max_steps << ")");
	    
	    if (verbose) {
		std::cout << "Minimize energy" << std::endl;
		std::cout << "  current force_max = " << sqrt(f_max_sqr_) << std::endl;
		std::cout << "  f_tol = " << sqrt(ftol_sqr_) << std::endl;
		std::cout << "  max_steps = " << max_steps << std::endl;
		std::cout << "  max displacement = " << max_displacement_ << std::endl;
	    }
	    int iters = 0;
	    for (; iters < max_steps && f_max_sqr_ > ftol_sqr_; iters++) {
		steepestDescentStep();
		
		dp_MAX += sqrt(dp_sqr_max_);
		
		resort_flag_ = dp_MAX > skin_half;
                LOG4ESPP_INFO(theLogger, "maxDist = " << dp_MAX << ", skin/2 = " << skin_half);
		
		if (resort_flag_) {
                    LOG4ESPP_INFO(theLogger, "Particles will be decomposed.");
		    dp_MAX = 0.;
		    storage.decompose();
                    LOG4ESPP_INFO(theLogger, "Particles have been decomposed.");
		    resort_flag_ = false;
		}
		
		updateForces();
		
		if (verbose)
		    std::cout << nstep_ << ": f_max^2=" << f_max_sqr_ << " max_dp^2=" << dp_sqr_max_ << std::endl;
		
		nstep_++;
	    }
	    
	    if (verbose) {
		std::cout << "Minimize energy finished" << std::endl;
		std::cout << "  current force_max = " << sqrt(f_max_sqr_) << std::endl;
		std::cout << "  run for steps = " << iters << std::endl;
		std::cout << "  max displacement^2 = " << dp_sqr_max_ << std::endl;
		if (f_max_sqr_ > ftol_sqr_) {
		    std::cout << "WARNING: the current max force is greater than the ftol=" << sqrt(ftol_sqr_);
		    std::cout << " The system might required additional run of energy minimization." << std::endl;
		}
	    }
	    retval = (f_max_sqr_ < ftol_sqr_);
	    
	    LOG4ESPP_INFO(theLogger,
			  "finished run, f_max_sqr_^2=" << f_max_sqr_ << " max_displ^2=" << dp_sqr_max_);
	    return retval;
	}
	
	
	void MinimizeEnergy::updateForces() {
	    LOG4ESPP_INFO(theLogger,
			  "update ghosts, calculate forces and collect ghost forces");
	    System& system = getSystemRef();
	    
	    system.storage->updateGhosts();
	    
	    LOG4ESPP_INFO(theLogger, "calculate forces");
	    // First initialize
	    CellList localCells = system.storage->getLocalCells();
	    for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
		cit->force() = 0.0;
	    }
	    
	    // Calculate force component from all interactions.
	    const InteractionList& srIL = system.shortRangeInteractions;
	    
	    for (size_t i = 0; i < srIL.size(); i++) {
		LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
		srIL[i]->addForces();
	    }
	    // Collect forces from ghost particles.
	    system.storage->collectGhostForces();
	    
	    // Get max force in the system.
	    LOG4ESPP_DEBUG(theLogger, "get max force in the system");
	    real f_max = -std::numeric_limits<real>::max();
	    CellList realCells = system.storage->getRealCells();
	    for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
		f_max = std::max(f_max, cit->force().sqr());
	    }
	    mpi::all_reduce(*system.comm, f_max, f_max_sqr_, boost::mpi::maximum<real>());
	}
	
	template <typename T> int sgn(T val) {
	    return (T(0) < val) - (val < T(0));
	}
	
	void MinimizeEnergy::steepestDescentStep() {
	    LOG4ESPP_INFO(theLogger, "steepestDescent single step");
	    System& system = getSystemRef();
	    
	    real f_sqr, dp, dp_sqr;
	    real f_max = sqrt(f_max_sqr_);
	    real dp_sqr_max = std::numeric_limits<real>::min();
	    
	    // Iterate over only real particles.
	    CellList realCells = system.storage->getRealCells();
	    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {

		if (variable_step_flag_) {
		    dp_sqr = 0.;
		    for (int i = 0; i < 3; i++) {   // Perhaps it can be done better.
			dp = max_displacement_ * cit->force()[i]/f_max;
			dp_sqr += dp*dp;
			
			// Update position component by dp.
			cit->position()[i] += dp;
		    }
		    dp_sqr_max = std::max(dp_sqr_max, dp_sqr);
		    
		} else {
		    dp_sqr = 0.;
		    for (int i = 0; i < 3; i++) {   // Perhaps it can be done better.
			dp = gamma_ * cit->force()[i];
			if (fabs(dp) > max_displacement_)
			    dp = sgn<real>(dp)*max_displacement_;
			dp_sqr += dp*dp;
			
			// Update position component by dp.
			cit->position()[i] += dp;
		    }
		    dp_sqr_max = std::max(dp_sqr_max, dp_sqr);
		}
	    }

	    LOG4ESPP_INFO(theLogger, "steepestDescentStep calculating dp_sqr_max");
	    mpi::all_reduce(*system.comm, dp_sqr_max, dp_sqr_max_, boost::mpi::maximum<real>());
	}
	
	void MinimizeEnergy::registerPython() {
	    using namespace espressopp::python;
	    
	    // Note: use noncopyable and no_init for abstract classes
	    class_<MinimizeEnergy, boost::noncopyable>
		("integrator_MinimizeEnergy", init<shared_ptr<System>, real, real, real, bool>())
		.add_property("f_max", &MinimizeEnergy::getFMax)
		.add_property("displacement", &MinimizeEnergy::getDpMax)
		.add_property("step", make_getter(&MinimizeEnergy::nstep_), make_setter(&MinimizeEnergy::nstep_))
		.def("run", &MinimizeEnergy::run);
	}
	
    }  // end namespace integrator
}  // end namespace espressopp
