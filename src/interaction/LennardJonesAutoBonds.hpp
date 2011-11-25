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

/*

      inline bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        bool found = false;
        Particle &pa;
        Particle &pb;
        int p1_count = 0;
        int p2_count = 0;
        if (dist.sqr() <= cutoffSqr) {
          switch (max_crosslinks) {
            case 1 :
              for (FixedPairList::PairList::Iterator it(*bondlist); it.isValid(); ++it) {
            	if ( ( p1.id() == it->first->id() || p1.id() == it->second->id() ) ||
            	     ( p2.id() == it->first->id() || p2.id() == it->second->id() ) ) {
            	  found = true;
            	}
              if (!found) {
                bondlist->add(p1.id(), p2.id());
              }
              break;
              };
            case 2 :
              if (p1.ghost() && p2.ghost()) {
            	  printf("ERROR: both particles are ghosts in _computeForce !!!\n");
            	  printf("       This should never happen !\n");
              }
              else
              if (p1.ghost()) {
            	  pa = p2;
            	  pb = p1;
              }
              else
              if (p2.ghost()) {
            	  pa = p1;
            	  pb = p2;
              }
              else {
            	  pa = p1;
            	  pb = p2;
              }
              for (FixedPairList::PairList::Iterator it(*bondlist); it.isValid(); ++it) {
             	if ( pa.id() == it->first->id() || pa.id() == it->second->id() ) {
             		++p1_count;
             	}
             	if ( p2.id() == it->first->id() || p2.id() == it->second->id() ) {
             		++p2_count;
             	}
              }
              if (p1_count + p2_count <= 1) {
                  bondlist->add(pa.id(), pb.id());
              }
              if (p1_count + p2_count == 1) {
              }
              break;
            default :
              printf("WARNING: currently crosslinking will only work for pair and triple joints. No crosslinking done.\n");
              break;
          }
        }
        return _computeForce(force, dist);
      }
*/

/* this version will not run correctly with more than 1 CPU */
      inline bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        bool found = false;
        long pa;
        long pb;
        bool add_second = false;
        int p1_count = 0;
        int p2_count = 0;
        bool pa_isghost = false;
        bool pb_isghost = false;
        if (dist.sqr() <= cutoffSqr) {
          switch (max_crosslinks) {
            case 1 :
              for (FixedPairList::PairList::Iterator it(*bondlist); it.isValid(); ++it) {
            	if ( ( p1.id() == it->first->id() || p1.id() == it->second->id() ) ||
            	     ( p2.id() == it->first->id() || p2.id() == it->second->id() ) ) {
            	  found = true;
            	}
              if (!found) {
                bondlist->add(p1.id(), p2.id());
              }
              break;
              };
            case 2 :
              for (FixedPairList::PairList::Iterator it(*bondlist); it.isValid(); ++it) {
             	if ( p1.id() == it->first->id() || p1.id() == it->second->id() ) {
             		++p1_count;
             		if (p1.id() == it->first->id()) {
             			pa = p2.id();
             			pa_isghost = p2.ghost();
             			pb = it->second->id();
             			pb_isghost = it->second->ghost();
             		} else {
             			pa = p2.id();
             			pa_isghost = p2.ghost();
             			pb = it->first->id();
             			pb_isghost = it->first->ghost();
             		}
             	}
             	if ( p2.id() == it->first->id() || p2.id() == it->second->id() ) {
             		++p2_count;
             		if (p2.id() == it->first->id()) {
             			pa = p1.id();
             			pa_isghost = p1.ghost();
             			pb = it->second->id();
             			pb_isghost = it->second->ghost();
             		} else {
             			pa = p1.id();
             			pa_isghost = p1.ghost();
             			pb = it->first->id();
             			pb_isghost = it->first->ghost();
             		}
             	}
              }
              if (p1_count + p2_count <= 1) {
                  bondlist->add(p1.id(), p2.id());
              }
              if (p1_count + p2_count == 1) {
            	if (!(pa_isghost && pb_isghost)) {
                    bondlist->add(pa, pb);
            	} else {
            		printf("WARNING: didn't add second crosslink, because both particles were ghosts");
            	}
              }
              break;
            default :
              printf("WARNING: currently crosslinking will only work for pair and triple joints. No crosslinking done.\n");
              break;
          }
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
