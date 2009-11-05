#include <iostream>
#include <cmath>
#include <vector>

#include "Particle.hpp"
#include "PeriodicBC.hpp"

namespace espresso {
  typedef std::vector<Particle> Cell;

  /** InteractionBase provides loop templates to compute
      forces and energies of the various interactions. */

  template< typename Derived >
  class InteractionBase {
  public:
    // full square over two cells
    virtual real 
    computeCellEnergies(BoundaryConditions &bc, 
			Cell &cell1, Cell &cell2) {
      real e = 0.0;
      for (int i = 0, endi = cell1.size(); i < endi; i++)
	for (int j = 0, endj = cell2.size(); j < endj; j++) {
	  Particle &p1 = cell1[i];
	  Particle &p2 = cell2[j];
	  real dist[3];
	  real distSqr;
	  bc.getMinimumImageVector(dist, &distSqr, p1.r.p, p2.r.p);
	  e += computeEnergy(p1, p2, dist, distSqr);
	}
      return e;
    }

    // half square over a single cell
    virtual real 
    computeCellEnergies(BoundaryConditions &bc, 
			Cell &cell) {
      real e = 0.0;
      for (int i = 0; i < cell.size(); i++)
	for (int j = 0; j < i; j++) {
	  Particle &p1 = cell[i];
	  Particle &p2 = cell[j];
	  real dist[3];
	  real distSqr;
	  bc.getMinimumImageVector(dist, &distSqr, p1.r.p, p2.r.p);
	  e += computeEnergy(p1, p2, dist, distSqr);
	}
      return e;
    }

    // full square over two cells
    virtual void 
    computeCellForces(BoundaryConditions &bc,
                      Cell &cell1, Cell &cell2) {
      for (int i = 0, endi = cell1.size(); i < endi; i++)
        for (int j = 0, endj = cell2.size(); j < endj; j++) {
          Particle &p1 = cell1[i];
          Particle &p2 = cell2[j];
          real dist[3];
          real distSqr;
          real force[3] = {0.0, 0.0, 0.0};
          bc.getMinimumImageVector(dist, &distSqr, p1.r.p, p2.r.p);
          computeForce(p1, p2, dist, distSqr, force);

          for (int k = 0; k < 3; i++) {
            p1.f.f[k] +=  force[k];
            p2.f.f[k] += -force[k];
          }
        }
    }
   
    // half square over a single cell
    virtual void
    computeCellForces(BoundaryConditions &bc,
                      Cell &cell) {
      for (int i = 0, endi = cell.size(); i < endi; i++)
        for (int j = 0; j < i; j++) {
          Particle &p1 = cell[i];
          Particle &p2 = cell[j];
          real dist[3];
          real distSqr;
          real force[3] = {0.0, 0.0, 0.0};
          bc.getMinimumImageVector(dist, &distSqr, p1.r.p, p2.r.p);
          computeForce(p1, p2, dist, distSqr, force);

          for (int k = 0; k < 3; i++) {
            p1.f.f[k] +=  force[k];
            p2.f.f[k] += -force[k];
          }
        }
    }

    // need loops for Verlet lists
    // and bonded interactions

    real computeEnergy(Particle& p1, Particle& p2, 
		       const real dist[3], real distSqr) {
      Parameters &params = getParameters(p1.type, p2.type);
      if (distSqr < params.getCutoffSqr()) {
	return static_cast< Derived* >(this)->computeEnergy(p1, p2, params, distSqr);
      } else return 0.0;
    }

    void computeForce(Particle& p1, Particle& p2, 
		      const real dist[3], real distSqr, 
		      real force[3]) {
      Parameters &params = getParameters(p1.type, p2.type);
      if (distSqr < params.getCutoffSqr()) {
	static_cast< Derived* >(this)->computeForce(p1, p2, params, dist, force);
      } else { 
	force[0] = force[1] = force[2] = 0.0; 
      }
    }

  protected:
    Derived::Parameters&
    getParameters(int type1, int type2);

    Derived::Parameters&
    createParameters(int type1, int type2);
  
    class ParametersBase {
    private:
      real cutoff;
      real cutoffSqr;

    public:
      ParametersBase() { setCutoff(0.0); }
      void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff * cutoff; }
      real getCutoff() const { return cutoff; }
      real getCutoffSqr() const { return cutoffSqr; }
    };
  };


  /** This class provides routines to compute forces and energies
      of the Lennard-Jones potential.

      \f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
      \left( \frac{\sigma}{r} \right)^{6} \right]
      \f]
  */

  class LennardJones : public InteractionBase< LennardJones > {
  public:
    LennardJones() {}

    void setParameters(int type1, int type2, real ep, real sg, real rc) {
      Parameters& params = createParameters(type1, type2);
      params.setEpsilon(ep);
      params.setSigma(sg);
      params.setCutoff(rc);
    }

  private:
    friend class InteractionBase<>;
    /* nested class to store coefficients */
    class Parameters : public ParametersBase {
      real epsilon;
      real sigma;

    public:
      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }
    };

    real computeEnergy(Particle &p1, Particle &p2, 
		       Parameters &params, const real dist[3], 
		       const real distSqr) const;

    void computeForce(Particle& p1, Particle& p2, 
		      Parameters &params, const real dist[3],
		      const real distSqr,
		      real force[3]) const;
  };

  real 
  LennardJones::
  computeEnergy(Particle &p1, Particle &p2, 
		Parameters &params, const real dist[3], 
		const real distSqr) const {
    real frac2 = params.getSigma() * params.getSigma() / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * params.getEpsilon() * (frac6 * frac6 - frac6);
    return energy;
  }

  void 
  LennardJones::
  computeForce(Particle& p1, Particle& p2, 
	       Parameters &params, const real dist[3],
	       const real distSqr, real force[3]) const {
    real frac2;
    real frac6;
    real distSqrInv;
    real ffactor;

    distSqrInv = 1.0 / distSqr;
    frac2 = params.getSigma() * params.getSigma() * distSqrInv;
    frac6 = frac2 * frac2 * frac2;
    ffactor = 48.0 * params.getEpsilon() * (frac6 * frac6 - 0.5 * frac6) * distSqrInv;

    for (int i = 0; i < 3; i++)
      force[i] = dist[i] * ffactor;
  }
}

int main() {

  /* create two particles and initialize their positions and forces */
  Particle p1;
  Particle p2;

  ParticlePosition r1 = {0.0, 0.0, 0.0};
  ParticlePosition r2 = {0.8, 1.1, 0.4};

  p1.r = r1;
  p2.r = r2;

  ParticleForce f1 = {0.0, 0.0, 0.0};
  ParticleForce f2 = {0.0, 0.0, 0.0};

  p1.f = f1;
  p2.f = f2;

  /* compute separation distance */
  Real3D dist = {p1.r.p[0] - p2.r.p[0], p1.r.p[1] - p2.r.p[1], p1.r.p[2] - p2.r.p[2]};
  real distSqr = dist.sqr();

  /* compute energy and force */
  InteractionBase<LennardJones> ib;
  ib.setParameters(1.0, 1.0, 2.0);
  std::cout << ib.computeEnergySqr(p1, p2, distSqr) << std::endl;
  ib.computeForce(p1, p2, dist);

  return 0;
}
