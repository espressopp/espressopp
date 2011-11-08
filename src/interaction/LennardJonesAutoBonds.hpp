// ESPP_CLASS
#ifndef _INTERACTION_LENNARDJONESAUTOBONDS_HPP
#define _INTERACTION_LENNARDJONESAUTOBONDS_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairList.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */
    class LennardJonesAutoBonds : public PotentialTemplate< LennardJonesAutoBonds > {
    private:
      real epsilon;
      real sigma;
      real ff1, ff2;
      real ef1, ef2;
      shared_ptr < FixedPairList > bondlist;

    public:
      static void registerPython();

      LennardJonesAutoBonds()
	: epsilon(0.0), sigma(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      LennardJonesAutoBonds(real _epsilon, real _sigma,
		   real _cutoff, real _shift, shared_ptr < FixedPairList > _fpl)
	: epsilon(_epsilon), sigma(_sigma), bondlist(_fpl) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      LennardJonesAutoBonds(real _epsilon, real _sigma,
		   real _cutoff, shared_ptr < FixedPairList > _fpl)
	: epsilon(_epsilon), sigma(_sigma), bondlist(_fpl) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
        setAutoShift(); 
      }

      void preset() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        updateAutoShift();
        preset();
      }
      
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { 
        sigma = _sigma; 
        updateAutoShift();
        preset();
      }
      real getSigma() const { return sigma; }

      real computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        return computeEnergy(dist);
      }

      void setCutoff(real _cutoff) {
        cutoff = _cutoff;
        cutoffSqr = cutoff*cutoff;
        updateAutoShift();
      }

      real setAutoShift() {
        autoShift = true;
        if (cutoffSqr == infinity)
  	    shift = 0.0;
        else
  	      shift = derived_this()->_computeEnergySqrRaw(cutoffSqr);
        return shift;
      }

      void updateAutoShift() {
         if (autoShift) setAutoShift();
      }

      real computeEnergy(const Real3D& dist) const {
  	     return computeEnergySqr(dist.sqr());
      }

      real computeEnergy(real dist) const {
         return computeEnergySqr(dist*dist);
      }

      real computeEnergySqr(real distsqr) const {
         return _computeEnergySqr(distsqr);
      }

      real _computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        return _computeEnergy(dist);
      }

      real _computeEnergy(const Real3D& dist) const {
        return _computeEnergySqr(dist.sqr());
      }

      real _computeEnergy(real dist) const {
        return _computeEnergySqr(dist*dist);
      }

      real _computeEnergySqr(real distSqr) const {
        if (distSqr > cutoffSqr)
          return 0.0;
        else {
          real e = derived_this()->_computeEnergySqrRaw(distSqr) - shift;
          LOG4ESPP_TRACE(theLogger, "Epot(r*r=" << distSqr << ") = " << e);
          return e;
        }
      }


      Real3D computeForce(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        return computeForce(dist);
      }

      Real3D computeForce(const Real3D& dist) const {
        Real3D force;
        if(!_computeForce(force, dist))
          force = 0.0;
        return force;
      }

      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        bool found = false;
        if (dist.sqr() <= cutoffSqr) {
          for (FixedPairList::PairList::Iterator it(*bondlist); it.isValid(); ++it) {
            if ( (it->first->id() == p1.id() && it->second->id() == p2.id()) ||
                 (it->first->id() == p2.id() && it->second->id() == p1.id()) ) {
              found = true;
            }
          }
          if ( !found ) {
            if (bondlist->add(p1.id(), p2.id())) {
              // printf("added particle pair (%i, %i) to bond list\n", p1.id(), p2.id());
              ;
            }
          }
        }
        return _computeForce(force, dist);
      }

      bool _computeForce(Real3D& force, const Real3D& dist) const {
        real distSqr = dist.sqr();
        if (distSqr > cutoffSqr)
          return false;
        else {
          return derived_this()->_computeForceRaw(force, dist, distSqr);
        }
      }

      real _computeEnergySqrRaw(real distSqr) const {
        real frac2 = sigma*sigma / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {

        real frac2 = 1.0 / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
        force = dist * ffactor;
        return true;
      }
    };
  }
}

#endif
