#include "BC.hpp"

using namespace espresso;
using namespace espresso::bc;

void 
BC::fold(const real &pos[3]) const {
  real res[3]= {pos[0], pos[1], pos[2]};
  foldThis(res);
}
