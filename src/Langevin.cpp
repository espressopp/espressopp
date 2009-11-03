#include <cmath>
#include <cstdio>

#define PMASS(pt) 1

#define COORD_FIXED(coord) (2L << coord)

typedef struct {
  double f[3];
} ParticleForce;

typedef struct {
  double v[3];
} ParticleMomentum;

typedef struct {
  int ext_flag;
} ParticleLocal; 

typedef struct {
  ParticleForce f;
  ParticleMomentum m;
  ParticleLocal l;
} Particle;

double d_random(void)
{
  return 1.0;
}

typedef struct {
  /** The particles payload */
  Particle *part;
  /** Number of particles contained */
  int n;
  /** Number of particles that fit in until a resize is needed */
  int max;
} ParticleList;

typedef ParticleList Cell;

typedef struct {
  Cell **cell;
  int n;
  int max;
} CellPList;

struct System {
  double timeStep;
  CellPList localCells;
};

/** The Langevin thermostat consists of a friction and noise term coupled
    via the fluctuation-dissipation theorem. The friction term is a
    function of the particle velocities. 

    Open question: If the feature \feature{ROTATION} is compiled in, the rotational
    degrees of freedom are also coupled to the thermostat.

*/

class Langevin {

 public:
  
  /** Constructor of a Langevin thermostat.

      \param _system System object to which thermostat belongs
      \param gamma Friction constant for the Langevin thermostat.
  */

  Langevin(double _gamma, double _temperature)

  {
    temperature = _temperature;
    gamma       = _gamma;
  }

  void connect(MDIntegrator &integrator) {
    integrator.startIntegration.connect(init);
    // APPLY?
    // heatUp?
    // coolDown?
  }

  /** Prefactors have to be calculated only once before the first iteration */

  void init(MDIntegrator &integrator)

  { // calculate the prefactors

    real timestep = integrator.getTimeStep();
    pref1 = -gamma * timestep;
    pref2 = sqrt(24.0 * temperature * gamma / timestep);
  }

  /** very nasty: if we recalculate force when leaving/reentering the integrator,
      a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
      numbers are drawn twice, resulting in a different variance of the random force.
      This is corrected by additional heat when restarting the integrator here.
      Currently only works for the Langevin thermostat, although probably also others
      are affected.
  */

  void heatUp()
 
  {
    pref2buffer = pref2;
    pref2       *= sqrt(3.0);
  }

  /** Opposite to heatUp */

  void coolDown()
  {
    pref2 = pref2buffer;
  }

  /** compute force for single particle, called for intialization of forces */

  void frictionThermo(Particle *p)

  {
     double massf = sqrt(PMASS(*p));

     for (int j = 0 ; j < 3 ; j++) {
    
       if (!(p->l.ext_flag & COORD_FIXED(j))) {
         p->f.f[j] += pref1*p->m.v[j]*PMASS(*p) + pref2*(d_random()-0.5)*massf;
       } 
     }
  }

  /** This routine will be called by the Verlet integrator */

  void apply(System &system) 

  { /** loop over all local particles */

    printf("start apply\n");

    CellPList& localCells = system.getStorage().getActiveCells();

    for (int c = 0; c < localCells.n; c++) {
      printf ("local cell %d\n", c);
      Cell* cell = localCells.cell[c];
      Particle* p  = cell->part;
      int np = cell->n;
      for (int i = 0; i < np; i++) {
        frictionThermo(&p[i]);
      }
    }
  }

 private:

  System* system;  // system to which the thermostat belongs

  double gamma;
  double temperature;

  double pref1;
  double pref2;
  double pref2buffer;

};

int main(int argc, char **argv)

{
  ParticleList myParticleList;

  myParticleList.max = 2;
  myParticleList.n   = 1;
  myParticleList.part = new Particle [myParticleList.max];

  Particle& myParticle = myParticleList.part[0];

  myParticle.m.v[0] = 0.5;
  myParticle.m.v[1] = 0.3;
  myParticle.m.v[2] = 0.1;
  myParticle.f.f[0] = 0.0;
  myParticle.f.f[1] = 0.0;
  myParticle.f.f[2] = 0.0;
 
  CellPList myCellList;

  myCellList.max = 3;
  myCellList.n   = 1;
  myCellList.cell = new ParticleList*[myCellList.max];
  myCellList.cell[0] = &myParticleList;

  System mySystem;

  mySystem.timeStep = 0.005;
  mySystem.localCells = myCellList;

  double gamma = 0.3;
  double temperature = 1.5;

  Langevin* thermostat = new Langevin(&mySystem, 0.3, 0.2);

  double *f = myParticle.f.f;

  printf("my old force = %f %f %f\n", f[0], f[1], f[2]);

  thermostat->init();

  for (int i = 0; i < 10; i++) {
     thermostat->apply();
  }

  printf("my new force = %f %f %f\n", f[0], f[1], f[2]);

}
