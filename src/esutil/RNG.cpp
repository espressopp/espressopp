#include "RNG.hpp"
#include "python.hpp"
#include "mpi.hpp"
#include "types.hpp"

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

shared_ptr< NormalVariate >
RNG::createNormalVariate(const real mean, const real sigma) {
  return
    make_shared< NormalVariate >(boostRNG, mean, sigma);
}

std::vector< real >
RNG::uniformOnSphere() {
  return uniformOnSphereVariate();
}

shared_ptr< UniformOnSphere >
RNG::createUniformOnSphere(int dim) {
  return
    make_shared< UniformOnSphere >(boostRNG, dim);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
RNG::registerPython() {
  using namespace espresso::python;

   real (RNG::*pyoperator)()=&RNG::operator();

   class_< RNG >("esutil_RNG")
    .def("__call__",pyoperator)
    .add_property("normal",&RNG::normal)
    ;
}
