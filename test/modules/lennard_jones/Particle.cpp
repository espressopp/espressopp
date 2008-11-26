//particle constructor and print function

#include <string>
#include "Particle.hpp"

Particle::Particle(double _x, double _y, double _z) {
  x = _x;
  y = _y;
  z = _z;
};

std::string Particle::toString() {
  std::string str;
  char buffer[50];
  int n;
  n = sprintf(buffer, "<%f, %f, %f>", x, y, z);
  str = buffer;
  return str;
};
