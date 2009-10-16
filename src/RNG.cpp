#include "RNG.hpp"

#include "mpi.hpp"

using namespace espresso;
using namespace boost;

RNG::NormalVariate::NormalVariate(shared_ptr< RNGType > _rng,
				  const real mean, 
				  const real sigma) 
  : Super(*_rng, DistType(mean, sigma)), rng(_rng)
{}

RNG::UniformOnSphere::UniformOnSphere(shared_ptr< RNGType > _rng,
				      int dim)
  : Super(*_rng, DistType(dim)), rng(_rng)
{}


RNG::RNG(long _seed) 
  : boostRNG(make_shared< RNGType >()), 
    normalVariate(boostRNG),
    uniformOnSphereVariate(boostRNG)
{
  seed(_seed);
}

// long
// RNG::seedNondet() {
//   // get a nondeterministic seed on the first proc, distribute it to
//   // the other procs, and seed
//   //   long ndSeed = nondet();
//   //   seed(ndSeed);
//   //   return ndSeed;
// }

void
RNG::seed(long _seed) {
  // Seed the RNG for the given CPU
  boostRNG->seed(_seed + mpiWorld.rank());
}

real
RNG::operator()() { 
  return uniform_01<>()(*boostRNG);
}

int
RNG::operator()(int N) { 
  return variate_generator< RNGType, uniform_smallint<> >
    (*boostRNG, uniform_smallint<>(0, N-1))();
}


real
RNG::normal() {
  return normalVariate();
}

RNG::NormalVariate::SelfPtr 
RNG::createNormalVariate(const real mean, const real sigma) {
  return
    make_shared< RNG::NormalVariate >(boostRNG, mean, sigma);
}

std::vector< real >
RNG::uniformOnSphere() {
  return uniformOnSphereVariate();
}

RNG::UniformOnSphere::SelfPtr
RNG::createUniformOnSphere(int dim) {
  return
    make_shared< RNG::UniformOnSphere >(boostRNG, dim);
}
