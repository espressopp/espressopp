#include "NormalVariate.hpp"

using namespace espresso::esutil;
using namespace boost;

NormalVariate::NormalVariate(shared_ptr< RNGType > _rng,
				  const real mean, 
				  const real sigma) 
  : Super(*_rng, DistType(mean, sigma)), rng(_rng)
{}
