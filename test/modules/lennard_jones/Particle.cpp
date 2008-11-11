//particle constructor and print function

#include <string>
#include "Particle.hpp"

Particle::Particle(double xx, double yy, double zz) {
  x = xx;
  y = yy;
  z = zz;
};

std::string Particle::toString() {
  std::string str;
  char buffer[50];
  int n;
  n = sprintf(buffer, "<%f, %f, %f>", x, y, z);
  str = buffer;
  return str;
};
