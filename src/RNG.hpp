#ifndef _RNG_HPP
#define _RNG_HPP
#include <boost/random.hpp>
#include <vector>

#include "types.hpp"

namespace espresso {
  // CONFIG START
  // This should be used, but currently, boost::random seems to be buggy
  typedef boost::lagged_fibonacci607 RNGType;

  // If you REALLY need speed!
  //typedef boost::rand48 RNGType;
  // CONFIG END

  /** A class that allows to choose between different RNGs, but that
  provides a uniform interface for these.
  */
  using namespace boost;

  class RNG {
    shared_ptr< RNGType > boostRNG;

  public:
    // NESTED CLASSES BEGIN
    /** Nested class to generates normal distributed random
	variates. The class also keeps a shared pointer to the boost
	RNG object, so that it is not destroyed if the espresso RNG
	object dies. */
    class NormalVariate 
      : variate_generator< RNGType, normal_distribution< real > > 
    {
      typedef normal_distribution< real > DistType;
      typedef variate_generator< RNGType, DistType > Super;

      /// store the shared pointer to the RNG
      shared_ptr< RNGType > rng;

    public:
      typedef espresso::shared_ptr< NormalVariate > SelfPtr;

      NormalVariate(shared_ptr< RNGType > _rng,
		    const real mean = 0.0, 
		    const real sigma = 1.0);
      using Super::operator();
    };

    /** Nested class to generate random vectors that are uniformly
	distributed on a sphere. */
    class UniformOnSphere
      : variate_generator< RNGType, uniform_on_sphere< real > > 
    {
      typedef uniform_on_sphere< real > DistType;
      typedef variate_generator< RNGType, DistType > Super;

      /// store the shared pointer to the RNG
      shared_ptr< RNGType > rng;

    public:
      typedef shared_ptr< UniformOnSphere > SelfPtr;

      UniformOnSphere(shared_ptr< RNGType > _rng, int dim = 3);
      using Super::operator();
    };
    // NESTED CLASSES END

  public:
    typedef shared_ptr< RNG > SelfPtr;

    /** Init the RNG, use the given seed. */
    RNG(long _seed = 12345);

//     /** Seed it with a non-deterministic seed. */
//     long seedNondet();

    /** Seed the RNG. */
    void seed(long _seed);
    
    /** returns a uniformly distributed random number between 0 and
	1. */
    real operator()();

    /** returns a uniformly distributed integer random number in the
	interval [0, N-1]. */
    int operator()(int N);

    /** returns a normal distributed random number, with mean of 0.0
	and sigma of 1.0. If you need to generate many normal
	idstributed random numbers with different mean or sigma,
	use createNormalVariate() to create a variate generator object. */
    real normal();

    /** returns a variate generator object that generates random
	variates with a normal distribution and the given parameters. */
    NormalVariate::SelfPtr 
    createNormalVariate(const real mean = 0.0,
			const real sigma = 1.0);

    /** returns a random 3D vector that is uniformly distributed on a sphere. */
    std::vector< real > uniformOnSphere();

    UniformOnSphere::SelfPtr
    createUniformOnSphere(int dim = 3);

  private:
    // the standard normal variate (mean 0, sigma 1)
    NormalVariate normalVariate;
    UniformOnSphere uniformOnSphereVariate;
  };
}

#endif
