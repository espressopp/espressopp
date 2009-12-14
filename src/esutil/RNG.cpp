#include "RNG.hpp"

#include "mpi.hpp"

using namespace espresso::esutil;
using namespace espresso;
using namespace boost;

RNG::RNG(long _seed) 
  : boostRNG(make_shared< RNGType >()), 
    normalVariate(boostRNG),
    uniformOnSphereVariate(boostRNG)
{
  seed(_seed);
}

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

NormalVariate::SelfPtr 
RNG::createNormalVariate(const real mean, const real sigma) {
  return
    make_shared< NormalVariate >(boostRNG, mean, sigma);
}

std::vector< real >
RNG::uniformOnSphere() {
  return uniformOnSphereVariate();
}

UniformOnSphere::SelfPtr
RNG::createUniformOnSphere(int dim) {
  return
    make_shared< UniformOnSphere >(boostRNG, dim);
}
