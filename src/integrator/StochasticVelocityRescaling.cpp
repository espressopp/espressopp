#include "python.hpp"
#include "StochasticVelocityRescaling.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include <math.h>

#define BOLTZMANN 1.0 //in reduced units
namespace espresso {

using namespace iterator;

namespace integrator {

LOG4ESPP_LOGGER(StochasticVelocityRescaling::theLogger,
		"StochasticVelocityRescaling");

StochasticVelocityRescaling::StochasticVelocityRescaling(
		shared_ptr<System> system) :
		SystemAccess(system) {
	temperature = 0.0;
	coupling = 1; //tau_t coupling

	if (!system->rng) {
		throw std::runtime_error("system has no RNG");
	}

	rng = system->rng;

	//chose on of currently three implementations of a gamma-distributed random number
	gammaDist = new GammaDistributionBoost(rng);
	gammaDist = new GammaDistributionNR2nd(rng);
	gammaDist = new GammaDistributionNR3rd(rng);
	//TODO benchmark those distributions! a first dirty benchmark did not show any significant difference in computational time

	LOG4ESPP_INFO(theLogger, "StochasticVelocityRescaling constructed");
}

void StochasticVelocityRescaling::setTemperature(real _temperature) {
	temperature = _temperature;
}

real StochasticVelocityRescaling::getTemperature() {
	return temperature;
}

void StochasticVelocityRescaling::setCoupling(real _coupling) {
	coupling = _coupling;
}

real StochasticVelocityRescaling::getCoupling() {
	return coupling;
}

StochasticVelocityRescaling::~StochasticVelocityRescaling() {
}

void StochasticVelocityRescaling::rescaleVelocities() {
	LOG4ESPP_DEBUG(theLogger, "rescaleVelocities");

	int NPart_local, NPart;
	real EKin = 0.0;
	real EKin_local = 0.0;
	real EKin_ref = 0.0;
	int DegreesOfFreedom;
	real EKin_new = 0.0;
	real ScalingFactor;

	System& system = getSystemRef();
	CellList realCells = system.storage->getRealCells();

	for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
		Real3D vel = cit->velocity();
		EKin_local += 0.5 * cit->mass() * (vel * vel); //FIXME do not forget that his is only correct for velocity-verlet; leap-frog would require adjustments
	}

	NPart_local = system.storage->getNRealParticles();

	boost::mpi::all_reduce(*getSystem()->comm, EKin_local, EKin,
			std::plus<real>());
	boost::mpi::all_reduce(*getSystem()->comm, NPart_local, NPart,
			std::plus<int>());

	DegreesOfFreedom = 3.0 * NPart; //TODO this is _only_ true for simple system without any constraints
	//calculate the reference kinetic energy based on reference temperature 'temperature'
	EKin_ref = 0.5 * temperature * BOLTZMANN * DegreesOfFreedom;
	EKin_new = stochasticVR_pullEkin(EKin, EKin_ref, DegreesOfFreedom, coupling,
			rng);
	// it should always be larger than 0
	if (EKin_new <= 0)
		throw std::runtime_error(
				"EKin_new in StochasticVelocityRescaling::rescaleVelocities() is equal or smaller than 0");

	ScalingFactor = sqrt(EKin_new / EKin);

	for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
		cit->velocity() *= ScalingFactor;
	}

}

const real StochasticVelocityRescaling::stochasticVR_sumGaussians(const int n) {

	/** plain implementation **/
	/*	real tmp, sum = 0.0;
	 for (int j = 0; j < n; j++) {
	 tmp = rng->normal();
	 sum += tmp * tmp;
	 }
	 return sum;*/

	real rr;
	if (n == 0)
		return 0.0;
	else if (n == 1) {
		rr = rng->normal();
		return rr * rr;
	} else if (n % 2 == 0) {
		return 2.0 * gammaDist->drawNumber(n / 2);
	} else {
		rr = rng->normal();
		GammaDistributionNR3rd
		return 2.0 * gammaDist->drawNumber((n - 1) / 2) + rr * rr;
	}
}

const real StochasticVelocityRescaling::stochasticVR_pullEkin(real Ekin,
		real Ekin_ref, int dof, real taut, shared_ptr<esutil::RNG> rng) {
	real factor, rr;

	/*
	 the parameter taut is the coupling time divided by the timestep, i.e.
	 even a strong coupling (low coupling time) with a large timestep should be
	 in the order of 0.1 / 0.01 = 10
	 make sure it's reasonably large
	 */
	if (taut < 0.1)
		throw std::runtime_error("taut in stochasticVR_pullEkin is very low."
				"If that is intended, change the code");
	else {
		factor = exp(-1.0 / taut);
	}
	rr = rng->normal();
	return Ekin
			+ (1.0 - factor)
					* (Ekin_ref * (stochasticVR_sumGaussians(dof - 1) + rr * rr)
							/ dof - Ekin)
			+ 2.0 * rr * sqrt(Ekin * Ekin_ref / dof * (1.0 - factor) * factor);
}

const real GammaDistributionBoost::drawNumber(const unsigned int ia) {
	return rng->gamma(ia);
}

/** Gamma distribution, from Numerical Recipes, 2nd edition, pages 292 & 293 */
const real GammaDistributionNR2nd::drawNumber(const unsigned int ia) {
	real x, v1, v2, y, am, s, e;

	if (ia < 1)
		throw std::runtime_error(
				"Error in routine stochasticVR_gammaDeviate2nd"); //TODO: how do I get __FUNCTION__ in here?
	if (ia < 6) {
		x = 1.0;
		for (int j; j <= ia; j++)
			x = x * (*rng)();
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1 = 2.0 * (*rng)() - 1.0;
					v2 = 2.0 * (*rng)() - 1.0;
				} while (v1 * v1 + v2 * v2 > 1.0);
				y = v2 / v1;
				am = ia - 1;
				s = sqrt(2.0 * am + 1.0);
				x = s * y + am;
			} while (x <= 0.0);
			e = (1.0 + y * y) * exp(am * log(x / am) - s * y);
		} while ((*rng)() > e);
	}
	return x;
}

/** Gamma distribution, from Numerical Recipes, 3rd edition, pages 370 & 371 */
const real GammaDistributionNR3rd::drawNumber(const unsigned int ia) {
	real a1, a2, x, v, u;
	if (ia <= 0)
		throw std::runtime_error(
				"Error in routine stochasticVR_gammaDeviate3rd");
	a1 = ia - 1.0 / 3.0;
	a2 = 1.0 / sqrt(9.0 * a1);
	do {
		do {
			x = rng->normal();
			v = 1.0 + a2 * x;
		} while (v <= 0.0);
		v = v * v * v;
		u = (*rng)();
	} while ((u > 1.0 - 0.331 * x * x * x * x)
			&& log(u) > 0.5 * x * x + a1 * (1.0 - v + log(v)));
	return a1 * v;
}

/****************************************************
 ** REGISTRATION WITH PYTHON
 ****************************************************/

void StochasticVelocityRescaling::registerPython() {

	using namespace espresso::python;

	class_<StochasticVelocityRescaling, shared_ptr<StochasticVelocityRescaling> >

	("integrator_StochasticVelocityRescaling", init<shared_ptr<System> >())

	.add_property("temperature", &StochasticVelocityRescaling::getTemperature,
			&StochasticVelocityRescaling::setTemperature).add_property(
			"coupling", &StochasticVelocityRescaling::getCoupling,
			&StochasticVelocityRescaling::setCoupling);
}

}
}

