#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "estypes.hpp"

struct System {
  real boxL[3];

  const real *getBoxL() const { return boxL; }
};

#endif
