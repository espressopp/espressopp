//method to compute the potential energy

#include <math.h>
#include "LennardJonesInteraction.hpp"

double LennardJonesInteraction::computeLennardJonesEnergy(double _rsq) {

  double frac2;
  double frac6;

  frac2 = 1.0 / _rsq;
  frac6 = frac2 * frac2 * frac2;
  return 4.0 * 1.0 * (pow(frac6, 2) - frac6);
};
