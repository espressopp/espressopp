//method to compute the potential energy

#include <math.h>
#include "LennardJonesInteraction.hpp"

double LennardJonesInteraction::computeLennardJonesEnergy(double dist2) {
  double frac2;
  double frac6;

  if(dist2 < rc2) {
    frac2 = 1.0 / dist2;
    frac6 = frac2 * frac2 * frac2;
    return 4.0 * 1.0 * (pow(frac6, 2) - frac6);
  }
};
