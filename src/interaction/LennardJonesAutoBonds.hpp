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

// ESPP_CLASS
#ifndef _INTERACTION_LENNARDJONESAUTOBONDS_HPP
#define _INTERACTION_LENNARDJONESAUTOBONDS_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairList.hpp"
#include "Potential.hpp"

namespace espressopp {
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
      int max_crosslinks;

    public:
      static void registerPython();

      LennardJonesAutoBonds()
	: epsilon(0.0), sigma(0.0) {
        setShift(0.0);
        setCutoff(0.0);
        max_crosslinks = 0;
      }

      LennardJonesAutoBonds(real _epsilon, real _sigma,
		   real _cutoff, real _shift, shared_ptr < FixedPairList > _fpl, int _max_crosslinks)
	: epsilon(_epsilon), sigma(_sigma), bondlist(_fpl), max_crosslinks(_max_crosslinks) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      LennardJonesAutoBonds(real _epsilon, real _sigma,
		   real _cutoff, shared_ptr < FixedPairList > _fpl, int _max_crosslinks)
	: epsilon(_epsilon), sigma(_sigma), bondlist(_fpl), max_crosslinks(_max_crosslinks) {
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

      real getMaxCrosslinks() const { return max_crosslinks; }

      void setMaxCrosslinks(int _max_crosslinks) {
    	  max_crosslinks = _max_crosslinks;
      }

      inline real computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        return computeEnergy(dist);
      }

      inline void setCutoff(real _cutoff) {
        cutoff = _cutoff;
        cutoffSqr = cutoff*cutoff;
        updateAutoShift();
      }

      inline real setAutoShift() {
        autoShift = true;
        if (cutoffSqr == infinity)
  	    shift = 0.0;
        else
  	      shift = derived_this()->_computeEnergySqrRaw(cutoffSqr);
        return shift;
      }

      inline void updateAutoShift() {
         if (autoShift) setAutoShift();
      }

      inline real computeEnergy(const Real3D& dist) const {
  	     return computeEnergySqr(dist.sqr());
      }

      inline real computeEnergy(real dist) const {
         return computeEnergySqr(dist*dist);
      }

      inline real computeEnergySqr(real distsqr) const {
         return _computeEnergySqr(distsqr);
      }

      inline real _computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        return _computeEnergy(dist);
      }

      inline real _computeEnergy(const Real3D& dist) const {
        return _computeEnergySqr(dist.sqr());
      }

      inline real _computeEnergy(real dist) const {
        return _computeEnergySqr(dist*dist);
      }

      inline real _computeEnergySqr(real distSqr) const {
        if (distSqr > cutoffSqr)
          return 0.0;
        else {
          real e = derived_this()->_computeEnergySqrRaw(distSqr) - shift;
          LOG4ESPP_TRACE(theLogger, "Epot(r*r=" << distSqr << ") = " << e);
          return e;
        }
      }

      inline Real3D computeForce(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        return computeForce(dist);
      }

      inline Real3D computeForce(const Real3D& dist) const {
        Real3D force;
        if(!_computeForce(force, dist))
          force = 0.0;
        return force;
      }

      inline bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
    	if (p1.ghost()) {
    		printf("WARNING (this should never happen !): didn't add crosslink, because particle1 is ghosts");
    		return _computeForce(force, dist);
    	}
        if (dist.sqr() <= cutoffSqr) {
          FixedPairList::GlobalPairs* globalPairs = bondlist->getGlobalPairs();
          if (globalPairs->count(p1.id()) + globalPairs->count(p2.id()) < max_crosslinks) {
       	    globalPairs->insert(globalPairs->begin(), std::make_pair(p1.id(), p2.id()));
            bondlist->add(p1.id(), p2.id());
          };
        }
        return _computeForce(force, dist);
      }

      inline bool _computeForce(Real3D& force, const Real3D& dist) const {
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
