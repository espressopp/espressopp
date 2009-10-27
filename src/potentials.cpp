#include <iostream>
#include <cmath>
#include <vector>

typedef double real;

/* from Particle.hpp */
struct ParticleProperties {
  size_t identity;
  size_t type;
};

struct ParticlePosition {
  real p[3];
};

struct ParticleForce {
  real f[3];
  ParticleForce operator* (real s) {
    ParticleForce tmp = {f[0] * s, f[1] * s, f[2] * s};
    return tmp;
  }
};

struct ParticleMomentum {
  real v[3];
};

struct ParticleLocal {
  size_t i[3];
  size_t ghost;
};

struct Particle {
  ParticleProperties p;
  ParticlePosition r;
  ParticleMomentum m;
  ParticleForce f;
  ParticleLocal l;
};

typedef std::vector<Particle> Cell;
/* above from Particle.hpp */

struct Real3D {
  real rij[3];
  real sqr() {
    return pow(rij[0], 2) + pow(rij[1], 2) + pow(rij[2], 2);
  }
  real operator [] (size_t i) { return rij[i]; }
};

/* CRTP */
template<typename Derived>
class InteractionBase {
public:

  /* the next method should be moved */
  void setParameters(real ep, real sg, real rc) {
    static_cast<Derived*>(this)->_setParameters(ep, sg, rc);
  }

  real computeEnergySqr(Particle p1, Particle p2, real distSqr) {
    static_cast<Derived*>(this)->_computeEnergySqr(p1, p2, distSqr);
  }

  void computeForce(Particle& p1, Particle& p2, Real3D dist) {
    static_cast<Derived*>(this)->_computeForce(p1, p2, dist);
  }

};

/* pairwise potential which implements computeEnergy and computeForce */
class LennardJones: public InteractionBase< LennardJones > {
private:

  /* nested class to store coefficients and cutoff */
  class Parameters {
  private:
    real epsilon;
    real sigma;
    real cutoff;
    real cutoffSqr;

  public:
    Parameters() {}
    Parameters(real ep, real sg, real rc) {
      setEpsilon(ep);
      setSigma(sg);
      setCutoff(rc);
    }

    void setEpsilon(real _epsilon) { epsilon = _epsilon; }
    real getEpsilon() const { return epsilon; }

    void setSigma(real _sigma) { sigma = _sigma; }
    real getSigma() const { return sigma; }

    void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff * cutoff; }
    real getCutoff() const { return cutoff; }

    real _getCutoffSqr() const { return cutoffSqr; }
  };

  Parameters params;

  public:

    LennardJones(real _epsilon, real _sigma, real _cutoff) {
      params = Parameters(_epsilon, _sigma, _cutoff);
    }

    void _setParameters(real ep, real sg, real rc) {
      params.setEpsilon(ep);
      params.setSigma(sg);
      params.setCutoff(rc);
    }

    real _computeEnergySqr(Particle p1, Particle p2, real distSqr);

    void _computeForce(Particle& p1, Particle& p2, Real3D dist) const;
};

real LennardJones::_computeEnergySqr(Particle p1, Particle p2, real distSqr) {
  if (distSqr < params._getCutoffSqr()) {
    real frac2 = params.getSigma() * params.getSigma() / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * params.getEpsilon() * (frac6 * frac6 - frac6);
    return energy;
  } else {
    return 0.0;
  }
}

void LennardJones::_computeForce(Particle& p1, Particle& p2, Real3D dist) const {
  real frac2;
  real frac6;
  real distSqrInv;
  real ffactor;
  real distSqr = dist.sqr();

  if (distSqr < params._getCutoffSqr()) {
    distSqrInv = 1.0 / distSqr;
    frac2 = params.getSigma() * params.getSigma() * distSqrInv;
    frac6 = frac2 * frac2 * frac2;
    ffactor = 48.0 * params.getEpsilon() * (frac6 * frac6 - 0.5 * frac6) * distSqrInv;

    ParticleForce f = {dist[0], dist[1], dist[2]};
    p1.f = f * ffactor; //need += operator
    p2.f = f * (-1.0 * ffactor);
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
