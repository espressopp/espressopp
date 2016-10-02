/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
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
#ifndef _INTEGRATOR_STOCHASTICVELOCITYRESCALING_HPP
#define _INTEGRATOR_STOCHASTICVELOCITYRESCALING_HPP

#include "types.hpp"
#include "logging.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"

namespace espressopp {
namespace integrator {

class GammaDistribution {
public:
	GammaDistribution(shared_ptr<esutil::RNG> _rng) :
			rng(_rng) {
	}
	virtual ~GammaDistribution() {
	}
	;
	virtual real drawNumber(const unsigned int alpha) = 0;

protected:
	shared_ptr<esutil::RNG> rng;
};

/** Gamma distribution, from Boost */
class GammaDistributionBoost: public GammaDistribution {
public:
	GammaDistributionBoost(shared_ptr<esutil::RNG> _rng) :
			GammaDistribution(_rng) {
	}
	real drawNumber(const unsigned int ia);
};

/** Gamma distribution, from Numerical Recipes, 2nd edition, pages 292 & 293 */
class GammaDistributionNR2nd: public GammaDistribution {
public:
	GammaDistributionNR2nd(shared_ptr<esutil::RNG> _rng) :
			GammaDistribution(_rng) {
	}
	real drawNumber(const unsigned int ia);
};

/** Gamma distribution, from Numerical Recipes, 3rd edition, pages 370 & 371 */
class GammaDistributionNR3rd: public GammaDistribution {
public:
	GammaDistributionNR3rd(shared_ptr<esutil::RNG> _rng) :
			GammaDistribution(_rng) {
	}
	real drawNumber(const unsigned int ia);
};

class StochasticVelocityRescaling: public Extension {

public:

	StochasticVelocityRescaling(shared_ptr<System> system);

	void setTemperature(real temperature);

	real getTemperature();

	void setCoupling(real coupling);

	real getCoupling();

	~StochasticVelocityRescaling();

	/** Sum n squared Gaussian numbers - shortcut via Gamma distribution */
	real stochasticVR_sumGaussians(const int n);

	/** Pull new value for the kinetic energy following the canonical distribution
	 *  Cite: Bussi et al JCP (2007) (there's a typo in the paper - this code is correct
	 *  Ekin: current kinetic energy
	 *  Ekin_ref: reference kinetic energy
	 *  dof: degrees of freedom
	 *  taut: coupling time/strength
	 *  */
	real stochasticVR_pullEkin(real Ekin, real Ekin_ref, int dof,
			real taut, shared_ptr<esutil::RNG> rng);

	/** Register this class so it can be used from Python. */
	static void registerPython();

private:
    boost::signals2::connection _runInit, _aftIntV;

	real temperature; //!< desired user temperature
	real coupling; // how strong is the coupling, i.e., tau_t coupling time
  real pref; // factor in the rescaling part

  int NPart, NPart_local, DegreesOfFreedom;
  real EKin_ref;

  void initialize();

	shared_ptr<esutil::RNG> rng; //!< random number generator

	GammaDistribution *gammaDist;

	void rescaleVelocities();

    void connect();
    void disconnect();

	/** Logger */
	static LOG4ESPP_DECL_LOGGER(theLogger);
};

}
}

#endif
