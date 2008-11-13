#define RAND_MAX
#define NUMBER_OF_PARTICLES 22000

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "Particle.hpp"
#include "LennardJonesInteraction.hpp"

int main() {

  //variables for storing random numbers
  double rx;
  double ry;
  double rz;

  //create a vector to store the particles
  std::vector<Particle> pc(NUMBER_OF_PARTICLES);

  //assign random positions to the particles on r[0, 20]
  for(int i = 0; i < pc.size(); i++) {
    rx = 20 * double(rand()) / RAND_MAX;
    ry = 20 * double(rand()) / RAND_MAX;
    rz = 20 * double(rand()) / RAND_MAX;
    pc[i] = Particle(rx, ry, rz);
  }

  //print a particle to standard output as a test
  std::cout << pc[1].toString() << std::endl;

  //create a LJ interaction and set its cutoff
  LennardJonesInteraction lj = LennardJonesInteraction();
  lj.setCutoff(2.5);

  //loop over all pairs assuming non-periodic BCs
  double en;
  double rsq;
  int j;

  en = 0.0;
  for(int i = 0; i < pc.size() - 1; i++) {
    Particle* Pi = &pc[i];
    for(j = i + 1; j < pc.size(); j++) {
      Particle* Pj = &pc[j];
      rsq   = pow(Pi->getx() - Pj->getx(), 2);
      rsq  += pow(Pi->gety() - Pj->gety(), 2);
      rsq  += pow(Pi->getz() - Pj->getz(), 2);
      en += lj.computeLennardJonesEnergy(rsq);
    }
  }

  //write out the total LJ energy
  std::cout << "en = " << en << std::endl;

  return 0;
}
