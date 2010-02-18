
#include "LennardJones.hpp"

#include "types.hpp"
#include <cassert>

using namespace espresso;
using namespace espresso::interaction;

/************************************************************************************/

LennardJones::LennardJones()
{
  ntypes = 0;
}

/************************************************************************************/

LennardJones::~LennardJones()
{
}

/************************************************************************************/

INTERACTION_ROUTINES(LennardJones)

/************************************************************************************/

void LennardJones::enableShift(int type1, int type2)
{
  Parameters& params = parameterArray[type1][type2];
  params.enableShift();
}

/************************************************************************************/

void LennardJones::setParameters(int type1, int type2, real epsilon, real sigma, real rc) 
{
  assert(type1 >= 0);
  assert(type2 >= 0);

  if ((type1 >= ntypes) || (type2 >= ntypes)) {

    // reallocate the type array

    int maxtypes = std::max(type1+1, type2+1);

    if (maxtypes > ntypes + 1) {
       //printf("WARNING: number of types increased from %d to %d\n",
       //                 ntypes, maxtypes);
    }

    if (maxtypes > 1) {
       //printf("ERROR: currently not supported\n");
       exit(-1);
    }
 
    // ToDo: reallocation

    ntypes = maxtypes;
  }

  Parameters& params = parameterArray[type1][type2];
  params.setEpsilon(epsilon);
  params.setSigma(sigma);
  params.setCutoff(rc);
  params.preset();
}

