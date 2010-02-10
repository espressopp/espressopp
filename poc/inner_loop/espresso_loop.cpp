#include <espresso_common.hpp>
#include <vector>
#include <iostream>

#include "Real3D.hpp"
#include "esutil/RNG.hpp"
#include "esutil/Timer.hpp"

using namespace espresso;
namespace espresso {
  const int NP0 = 10000;
  const int NP1 = 10000;
  const real L = 100.0;

  longint counter;

  struct Particle {
    Real3D pos;
    Real3D f;
    int type;
  };

  std::vector< Particle > particles;

  void setupSystem() {
    std::cout << "Setting up system..." << std::endl;
    esutil::RNG rng;
    particles.resize(NP0+NP1);
    // create particles
    int pid = 0;
    for (int i = 0; i < NP0; i++) {
      particles[pid].pos = Real3D(rng()*L, rng()*L, rng()*L);
      particles[pid].type = 0;
      ++pid;
    }
    for (int i = 0; i < NP1; i++) {
      particles[pid].pos = Real3D(rng()*L, rng()*L, rng()*L);
      particles[pid].type = 1;
      ++pid;
    }
    std::cout << "System set up." << std::endl;
  }

  //////////////////////////////////////////////////
  // ESPRESSO LOOP
  //////////////////////////////////////////////////
  struct LJParams {
    real epsilon;
    real sigma;
    real ff1, ff2;
    real ef1, ef2;

    LJParams() {}
  
    LJParams(const real _epsilon, const real _sigma) 
      : epsilon(_epsilon), sigma(_sigma)
    {
      real sig2 = sigma * sigma;
      real sig6 = sig2 * sig2 * sig2;
      ff1 = 48.0 * epsilon * sig6 * sig6;
      ff2 = 24.0 * epsilon * sig6;
      ef1 =  4.0 * epsilon * sig6 * sig6;
      ef2 =  4.0 * epsilon * sig6;
    }
  };
  
  LJParams ljParams[2][2];

  void addForce(Particle &p1, Particle &p2) {
    ++counter;
    LJParams &params = ljParams[p1.type][p2.type];
    Real3D f = p1.pos - p2.pos;
    real frac2 = 1.0 / f.sqr();
    real frac6   = frac2 * frac2 * frac2;
    real ffactor = frac6 * (params.ff1 * frac6 - params.ff2) * frac2;
    f *= ffactor;
    p1.f += f;
    p2.f -= f;
  }

  void espressoLoop() {
    std::cout << "Doing ESPResSo loop..." << std::endl;

    // create LJ parameters
    ljParams[0][0] = LJParams(1.0, 1.0);
    ljParams[1][0] = LJParams(1.5, 1.0);
    ljParams[0][1] = LJParams(1.5, 1.0);
    ljParams[1][1] = LJParams(2.0, 1.0);

    counter = 0;
    esutil::WallTimer timer;

    // now loop
    timer.reset();
    for (int i = 0; i < NP0+NP1; ++i) {
      Particle &p1 = particles[i];
      for (int j = 0; j < i; ++j) {
	Particle &p2 = particles[j];
	addForce(p1, p2);
      }
    }
    std::cout << "ESPResSo loop: " << timer << " counter=" << counter << std::endl;
  }


  //////////////////////////////////////////////////
  // INTERACTION-BASED LOOP
  //////////////////////////////////////////////////
  class PairComputer {
  public:
    virtual void compute(Particle &p1, Particle &p2) = 0;
  };

  class LennardJones : public PairComputer {
    real epsilon;
    real sigma;
    real ff1, ff2;
    real ef1, ef2;

  public:
    LennardJones() {}
    
    LennardJones(const real _epsilon, const real _sigma) 
      : epsilon(_epsilon), sigma(_sigma)
    {
      real sig2 = sigma * sigma;
      real sig6 = sig2 * sig2 * sig2;
      ff1 = 48.0 * epsilon * sig6 * sig6;
      ff2 = 24.0 * epsilon * sig6;
      ef1 =  4.0 * epsilon * sig6 * sig6;
      ef2 =  4.0 * epsilon * sig6;
    }

    virtual void compute(Particle &p1, Particle &p2) {
      ++counter;
      Real3D f = p1.pos - p2.pos;
      real frac2 = 1.0 / f.sqr();
      real frac6   = frac2 * frac2 * frac2;
      real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
      f *= ffactor;
      p1.f += f;
      p2.f -= f;
    }
  };

  typedef std::vector< Particle* > ParticleList;

  //////////////////////////////////////////////////
  // VERSION 1
  class PairSet {
  public:
    virtual void foreach(shared_ptr< PairComputer > computer) = 0;
  };

  class FullSquarePairSet : public PairSet {
    shared_ptr< ParticleList > p1;
    shared_ptr< ParticleList > p2;
    
  public:
    FullSquarePairSet(shared_ptr< ParticleList > _p1, 
		      shared_ptr< ParticleList > _p2)
      : p1(_p1), p2(_p2) 
    {}

    virtual void foreach(shared_ptr< PairComputer > computer) {
      ParticleList::iterator pit1 = p1->begin();
      ParticleList::iterator pend1 = p1->end();
      for (; pit1 != pend1; ++pit1) {
	ParticleList::iterator pit2 = p2->begin();
	ParticleList::iterator pend2 = p2->end();
	for (; pit2 != pend2; ++pit2)
	  computer->compute(**pit1, **pit2);
      }
    }
  };

  class HalfSquarePairSet : public PairSet {
    shared_ptr< ParticleList > p;

  public:
    HalfSquarePairSet(shared_ptr< ParticleList > _p)
      : p(_p)
    {}
    
    virtual void foreach(shared_ptr< PairComputer > computer) {
      ParticleList::iterator pbegin = p->begin();
      ParticleList::iterator pend = p->end();
      for (ParticleList::iterator pit1 = pbegin; pit1 != pend; ++pit1)
	for (ParticleList::iterator pit2 = pbegin; pit2 != pit1; ++pit2)
	  computer->compute(**pit1, **pit2);
    }
  };

  class Interaction {
    shared_ptr< PairSet > pairs;
    shared_ptr< PairComputer > computer;

  public:
    Interaction(shared_ptr< PairSet > _pairs, 
		shared_ptr< PairComputer > _computer) 
      : pairs(_pairs), computer(_computer)
    {}

    void addForces() {
      pairs->foreach(computer);
    }
  };

  typedef std::vector< Interaction > InteractionList;

  void interactionLoop() {
    std::cout << "Doing Interaction-based loop..." << std::endl;

    // create particle lists
    shared_ptr< ParticleList > type0 = make_shared< ParticleList >();
    shared_ptr< ParticleList > type1 = make_shared< ParticleList >();
    for (int pid = 0; pid < NP0+NP1; pid++)
      if (particles[pid].type == 0)
	type0->push_back(&particles[pid]);
      else
	type1->push_back(&particles[pid]);

    // create PairSets
    shared_ptr< HalfSquarePairSet > type0half = make_shared< HalfSquarePairSet >(type0);
    shared_ptr< HalfSquarePairSet > type1half = make_shared< HalfSquarePairSet >(type1);
    shared_ptr< FullSquarePairSet > type01full = make_shared< FullSquarePairSet >(type0, type1);

    // create computers
    shared_ptr< LennardJones > lj00 = make_shared< LennardJones >(1.0, 1.0);
    shared_ptr< LennardJones > lj10 = make_shared< LennardJones >(1.0, 1.5);
    shared_ptr< LennardJones > lj11 = make_shared< LennardJones >(1.0, 2.0);

    // create Interactions
    InteractionList interactionList;
    interactionList.push_back(Interaction(type0half, lj00));
    interactionList.push_back(Interaction(type1half, lj11));
    interactionList.push_back(Interaction(type01full, lj10));

    counter = 0;
    esutil::WallTimer timer;
    // now loop
    timer.reset();
    InteractionList::iterator it = interactionList.begin();
    InteractionList::iterator end = interactionList.end();
    for (; it != end; it++)
      it->addForces();
    std::cout << "Interaction-based loop: " << timer << " counter=" << counter << std::endl;
  }
}

int main() {
  initMPIEnv();

  espresso::setupSystem();
  espresso::espressoLoop();
  espresso::interactionLoop();
  
  finalizeMPIEnv();
  return 0;
}
