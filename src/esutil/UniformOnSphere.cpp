#include "UniformOnSphere.hpp"

using namespace espresso::esutil;
using namespace boost;

UniformOnSphere::UniformOnSphere(shared_ptr< RNGType > _rng,
				 int dim)
  : Super(*_rng, DistType(dim)), rng(_rng)
{}
