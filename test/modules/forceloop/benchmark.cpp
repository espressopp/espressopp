/** @file benchmark.cpp

    A simple benchmark of the all-pairs loop operating on an unsorted
    ParticleStorage.
 */

#define LOG4ESPP_LEVEL_WARN

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "espresso_common.hpp"
#include "types.hpp"
#include "storage/Storage.hpp"
#include "particles/Computer.hpp"
#include "pairs/All.hpp"
#include "potential/ForceComputer.hpp"
#include "bc/PeriodicBC.hpp"
#include "potential/LennardJones.hpp"
#include "esutil/Timer.hpp"

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::storage;
using namespace boost;
using namespace std;

/// We set up a system of N^3 particles on an NxNxN lattice in a cubic
/// box of size. The particles are set in a cubic box with side length
/// variance around the lattice grid points.

/// number of particles in each dimension
const int N = 20;
/// dimension of the cubic simulation box
const real size = 21.0;
/// how far the particles maximally deviate from the grid point
const real variance = 0.1;

const size_t repeatParticleLoop = 10000;

size_t npart = 0;

static inline real dround(real x) { return floor(x + 0.5); }

/***********************************************************************************
 * Testing the Espresso way of calculating things
 ***********************************************************************************/

class TestEspresso {
public:
  Storage::SelfPtr storage;
  Property< Real3D >::SelfPtr position, force;
  
  TestEspresso() {
    storage = make_shared< Storage >();
    position = make_shared< Property< Real3D > >(storage);
    force = make_shared< Property< Real3D > >(storage);
  }

  void addParticles(vector< Real3D > &positions);

  void calculateForces(real epsilon, real sigma, real cutoff) NOINLINE;
  
  void runEmptyPairLoop() NOINLINE;
  
  real calculateMinDist() NOINLINE;
  
  void runEmptyLoop() NOINLINE;
  
  real calculateAverage() NOINLINE;
  
  Real3D getForce(size_t i) {
    return force->at(ParticleId(i));
  }

  Real3D getPosition(size_t i) {
    return position->at(ParticleId(i));
  }
};

void TestEspresso::addParticles(vector< Real3D > &positions) {
  size_t cnt = 0;
  BOOST_FOREACH(Real3D pos, positions)
    {
      ParticleId id = ParticleId(cnt);
      storage->addParticle(id);
      (*position)[id] = pos;
      (*force)[id] = 0.0;
      cnt++;
    }
}

void TestEspresso::calculateForces(real epsilon, real sigma, real cutoff) {
  bc::PeriodicBC::SelfPtr pbc = make_shared< bc::PeriodicBC >(size);
  pairs::All::SelfPtr allpairs = make_shared< pairs::All >(pbc, storage, position);
  potential::LennardJones::SelfPtr ljint 
    = make_shared< potential::LennardJones >(epsilon, sigma, cutoff);
  pairs::Computer::SelfPtr forceCompute 
    = ljint->createForceComputer(force);
  allpairs->foreach(forceCompute);
}

class EmptyPairComputer: public pairs::Computer {
public:
  typedef shared_ptr< EmptyPairComputer > SelfPtr;
  EmptyPairComputer() {}

  void prepare(Storage::SelfPtr storage1,
	       Storage::SelfPtr storage2) {}

  virtual void apply(const Real3D dist,
		     ParticleHandle ref1,
		     ParticleHandle ref2) {}
};

void TestEspresso::runEmptyPairLoop() {
  bc::PeriodicBC::SelfPtr pbc = make_shared< bc::PeriodicBC >(size);
  pairs::All::SelfPtr allpairs = make_shared< pairs::All >(pbc, storage, position);
  pairs::Computer::SelfPtr ljc = make_shared< EmptyPairComputer >();
  allpairs->foreach(ljc);
}

class MinDistComputer: public pairs::Computer {
public:
  typedef shared_ptr< MinDistComputer > SelfPtr;
  real min;
  
  MinDistComputer(): min(1.0e10) {}
  
  void prepare(Storage::SelfPtr storage1,
	       Storage::SelfPtr storage2) {}

  virtual void apply(const Real3D dist,
		     ParticleHandle ref1,
		     ParticleHandle ref2) {
    real d = dist.sqr();
    if (min > d) {
      min = d;
    }
  }
};

real TestEspresso::calculateMinDist() {
  bc::PeriodicBC::SelfPtr pbc = make_shared< bc::PeriodicBC >(size);
  pairs::All::SelfPtr allpairs = make_shared< pairs::All >(pbc, storage, position);
  MinDistComputer::SelfPtr mincomp = make_shared< MinDistComputer >();
  allpairs->foreach(*mincomp);
  return sqrt(mincomp->min);
}

class AverageComputer: public particles::Computer {
public:
  typedef shared_ptr< AverageComputer > SelfPtr;

  PropertyHandle< Real3D > property;
  real average;

  AverageComputer(const PropertyHandle< Real3D > &_property)
    : property(_property), average(0) {}

  void prepare(Storage::SelfPtr storage) {}

  void apply(ParticleHandle ref) {
    const Real3D p = property[ref];
    average += sqrt(p.sqr());
  }
};

real TestEspresso::calculateAverage() {
  AverageComputer::SelfPtr avgcompute 
    = make_shared< AverageComputer >(force->getHandle(storage));
  storage->foreach(*avgcompute);
  return avgcompute->average;
}

class EmptyComputer: public particles::Computer {
public:
  EmptyComputer() {}
  void prepare(Storage::SelfPtr storage) {}
  virtual void apply(ParticleHandle ref) {}
};

void TestEspresso::runEmptyLoop() {
  EmptyComputer compute;
  storage->foreach(compute);
}

/***********************************************************************************
 * Testing the hardcoded way
 ***********************************************************************************/

class TestBasic {
public:
  vector<Real3D> position, force;
  
  TestBasic() {}

  void addParticles(vector< Real3D > &positions);
  
  void calculateForces(real epsilon, real sigma, real cutoff) NOINLINE;
  
  real calculateMinDist() NOINLINE;
  
  real calculateAverage() NOINLINE;
  
  Real3D getForce(size_t i) {
    return force[i];
  }
  
  Real3D getPosition(size_t i) {
    return position[i];
  }
};

void TestBasic::addParticles(vector< Real3D > &positions) {
  BOOST_FOREACH(Real3D pos, positions)
    {
      position.push_back(pos);
      force.push_back(0.0);
    }
}

void TestBasic::calculateForces(real epsilon, real sigma, real cutoff) {
  real cutoffSqr = cutoff*cutoff;
  real sizeInverse = 1./size;

  for (size_t i = 0; i < npart; ++i) {
    for (size_t j = i+1; j < npart; ++j) {
      Real3D pos1 = position[i];
      Real3D pos2 = position[j];
      Real3D dist = pos1 - pos2;

      for(size_t c = 0; c < 3; ++c) {
	dist[c] -= dround(dist[c]*sizeInverse)*size;
      }

      {
	real frac2;
	real frac6;
	real distSqrInv;
	real distSqr = dist.sqr();
                
	if (distSqr < cutoffSqr) {
	  Real3D f(0.0, 0.0, 0.0);
	  distSqrInv = 1.0 / distSqr;
	  frac2 = sigma*sigma * distSqrInv;
	  frac6 = frac2 * frac2 * frac2;
	  real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5*frac6) * distSqrInv;
	  f = dist * ffactor;
	  force[i] += f;
	  force[j] -= f;
	} 
      }
    }
  }
}

real TestBasic::calculateMinDist() {
    real min = 1e10;
    real sizeInverse = 1./size;

    for (size_t i = 0; i < npart; ++i) {
        for (size_t j = i+1; j < npart; ++j) {
            Real3D pos1 = position[i];
            Real3D pos2 = position[j];

            Real3D dist = pos1 - pos2;
            for(size_t c = 0; c < 3; ++c) {
                dist[c] -= dround(dist[c]*sizeInverse)*size;
            }
            real d = dist.sqr();

	    if (min > d) {
		min = d;
            }
        }
    }
    return sqrt(min);
}

real TestBasic::calculateAverage() {
    real average = 0;
    for (size_t i = 0; i < npart; ++i) {
        Real3D p = force[i];
        average += sqrt(p.sqr());
    }
    return average;
}

/***********************************************************************************
 * Testing environment
 ***********************************************************************************/

void generateParticles(vector< Real3D >& positions) {
  srand48(123);
  real grid = size / (N+1);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) { 
      for (int k = 0; k < N; k++) {
	Real3D pos(i*grid + drand48()*variance,
		   j*grid + drand48()*variance,
		   k*grid + drand48()*variance);
	positions.push_back(pos);
      }
    }
  }
//   positions.push_back(Real3D(0.0, 0.0, 0.0));
//   positions.push_back(Real3D(1.0, 0.0, 0.0));
}

/// this routine runs the tests as defined above
int main()
{
  LOG4ESPP_CONFIGURE();

  cout << "BENCHMARK" << endl;
  cout << "Computing LJ forces in a cubic box (size=" << size << ") of " 
       << N << "x" 
       << N << "x" 
       << N << " particles." << endl;
  cout << endl;

  esutil::WallTimer timer;
  TestEspresso espresso;
  cout << setw(30) << "setup Espresso: " << timer << endl;

  timer.reset();
  TestBasic basic;
  cout << setw(30) << "setup Basic: " << timer << endl;

  // generate particles in the particle storage

  vector< Real3D > positions;
  generateParticles(positions);
  npart = positions.size();

  timer.reset();
  espresso.addParticles(positions);
  cout << setw(30) << "generate Espresso: " << timer << endl;

  timer.reset();
  basic.addParticles(positions);
  cout << setw(30) << "generate Basic: " << timer << endl;

  ////////////////////////////////////////////////////////////////////////
  cout << endl << "PARTICLE PAIR LOOPING TESTS" << endl;
  ////////////////////////////////////////////////////////////////////////

  // EMPTY PAIR LOOP
  timer.reset();
  espresso.runEmptyPairLoop();
  cout << setw(30) << "empty pair loop: " << timer << endl << endl;

  // CALCULATE MINIMUM DISTANCE
  timer.reset();
  real minb = basic.calculateMinDist();
  real basicMinTime = timer.getElapsedTime();
  cout << setw(30) << "min Basic: " << timer << " (min=" << minb << ")" << endl; 

  timer.reset();
  real mine = espresso.calculateMinDist();
  real espressoMinTime = timer.getElapsedTime();
  cout << setw(30) << "min Espresso: " << timer << " (min=" << mine << ")" << endl; 

  cout << setw(30) << "RATIO: " << (espressoMinTime / basicMinTime) << endl;

  // check consistency
  if (abs(mine-minb)/abs(mine) > 1e-5)
    cout << "ERROR: minima are different: " << mine << " != " << minb << endl;
  cout << endl;

  // CALCULATE FORCES
  timer.reset();
  basic.calculateForces(1.1, 1.2, 2.5);
  real basicCalcTime = timer.getElapsedTime();
  cout << setw(30) << "calc basic: " << timer << endl;

  timer.reset();
  espresso.calculateForces(1.1, 1.2, 2.5);
  real espressoCalcTime = timer.getElapsedTime();
  cout << setw(30) << "calc Espresso: " << timer << endl;

  cout << setw(30) << "RATIO: " << (espressoCalcTime / basicCalcTime) << endl;

  // check consistency
  for (size_t i = 0; i < npart; ++i) {
    Real3D fe = espresso.getForce(i);
    Real3D fb = basic.getForce(i);
    real diff = sqrt((fe-fb).sqr());
    if (diff/sqrt(fe.sqr()) > 1.0e-5) {
      Real3D pe = espresso.getPosition(i);
      Real3D pb = basic.getPosition(i);
      cout << "ERROR: difference " << diff << " too big for particle " << i << endl;
      cout << "     basic: f=(" << fb[0] << ", " << fb[1] << ", " << fb[2] << ")" << endl;
      cout << "            p=(" << pb[0] << ", " << pb[1] << ", " << pb[2] << ")" << endl;
      cout << "  espresso: f=(" << fe[0] << ", " << fe[1] << ", " << fe[2] << ")" << endl;
      cout << "            p=(" << pe[0] << ", " << pe[1] << ", " << pe[2] << ")" << endl;
      cout << endl;
    }
  }
  cout << endl;

  ////////////////////////////////////////////////////////////////////////
  cout << endl << "PARTICLE LOOPING TESTS" << endl;
  ////////////////////////////////////////////////////////////////////////

  // EMPTY LOOP
  timer.reset();
  for (size_t cnt = 0; cnt < repeatParticleLoop; ++cnt)
    espresso.runEmptyLoop();
  cout << setw(30) << "empty loop: " << timer << endl << endl;

  // CALCULATE AVERAGE
  timer.reset();
  real avb;
  for (size_t cnt = 0; cnt < repeatParticleLoop; ++cnt)
    avb = basic.calculateAverage();
  real basicAvgTime = timer.getElapsedTime();
  cout << setw(30) << "average Basic: " << timer << " (avg=" << avb << ")" << endl;

  timer.reset();
  real ave;
  for (size_t cnt = 0; cnt < repeatParticleLoop; ++cnt)
    ave = espresso.calculateAverage();
  real espressoAvgTime = timer.getElapsedTime();
  cout << setw(30) << "average Espresso: " << timer << " (avg=" << ave << ")" << endl;

  cout << setw(30) << "RATIO: " << (espressoAvgTime / basicAvgTime) << endl;

  if (abs(ave-avb)/abs(ave) > 1e-5)
    cout << "ERROR: averages are different: " << ave << " != " << avb << endl;
  cout << endl;
}
