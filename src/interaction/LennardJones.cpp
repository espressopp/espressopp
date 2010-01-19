
#include "LennardJones.hpp"

#include "types.hpp"

using namespace espresso;
using namespace espresso::interaction;

LennardJones::LennardJones()
{
}

LennardJones::~LennardJones()
{
}

INTERACTION_ROUTINES(LennardJones)

void LennardJones::setParameters(int type1, int type2, real ep, real sg, real rc) 
{
  Parameters& params = parameterArray[type1][type2];
  params.setEpsilon(ep);
  params.setSigma(sg);
  params.setCutoff(rc);
}

