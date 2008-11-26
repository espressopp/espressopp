#define RAND_MAX

//NUMBER_OF_PARTICLES MUST BE A PERFECT CUBE
#define NUMBER_OF_PARTICLES 27000

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "Particle.hpp"
#include "LennardJonesInteraction.hpp"
#include "FullPairIterator.hpp"
#include "BoundaryConditions.hpp"

template<class Interaction, class PairIterator>
void compute(std::vector<Particle>& _pc, Interaction lj, PairIterator it, Cube _bc) {
  
  double rsq;
  double en;

  en = 0.0;

  //loop over all pairs assuming non-periodic BCs
  for (it.reset(); !it.done(); it.next()) {
    Particle Pi = _pc[it.first()];
    Particle Pj = _pc[it.second()];
    rsq=_bc.getMinimumImageDistance(Pi, Pj);
    /*
    rsq  = pow(Pi.getx() - Pj.getx(), 2);
    rsq += pow(Pi.gety() - Pj.gety(), 2);
    rsq += pow(Pi.getz() - Pj.getz(), 2);*/
    if(rsq < lj.getCutoffSq())
      en  += lj.computeLennardJonesEnergy(rsq);
  }

  //write out the total LJ energy
  std::cout << "en = " << (en / NUMBER_OF_PARTICLES) << std::endl;
}

void init_configuration(std::vector<Particle>& _pc, double _rho, double _dr, Cube& _bc) {

  double one_third = 1.0 / 3.0;
  int n_cube_root = static_cast<int>(pow(NUMBER_OF_PARTICLES, one_third) + 1);

  if(static_cast<int>(pow(n_cube_root, 3) != NUMBER_OF_PARTICLES)) {
    std::cout << "Mismatch error in N." << std::endl;
  }

  //variables for storing random numbers
  double rx;
  double ry;
  double rz;

  //particles are initialized on a simple cube lattice
  //with a lattice spacing of 1
  int iParticle = 0;
  for(int i = 0; i < n_cube_root; i++) {
    rx = i + 0.5;
    for(int j = 0; j < n_cube_root; j++) {
      ry = j + 0.5;
      for(int k = 0; k < n_cube_root; k++) {
        rz = k + 0.5;
        _pc[iParticle] = Particle(rx, ry, rz);
        iParticle++;
      }
    }
  }

  double box_length = pow(NUMBER_OF_PARTICLES / _rho, one_third);
  double scale_factor = box_length / n_cube_root;
  std::cout << "L=" << box_length << std::endl;

  //set the size of the cubic box
  _bc.setSide(box_length);

  iParticle = 0;
  double rand_disp;
  for(int i = 0; i < NUMBER_OF_PARTICLES; i++) {
    rand_disp = _dr * (2.0 * rand() / RAND_MAX - 1.0);
    _pc[iParticle].setx(scale_factor * _pc[iParticle].getx() + rand_disp);
    _pc[iParticle].sety(scale_factor * _pc[iParticle].gety() + rand_disp);
    _pc[iParticle].setz(scale_factor * _pc[iParticle].getz() + rand_disp);
    iParticle++;
  }
}

int main() {

  //specify the fluid density in LJ reduced units
  double rho = 0.6;

  //create a vector to store the particles where
  //number of particles must be the cube of an integer
  std::vector<Particle> pc(NUMBER_OF_PARTICLES);

  //create an object to handle the periodic boundary conditions
  Cube bc = Cube();

  //initialize particle on a simple cubic lattice
  //with slight random displacements of max(-dr, dr)
  double dr = 0.1;
  init_configuration(pc, rho, dr, bc);

  //print the first, second and last particle coordinates
  std::cout << pc[0].toString() << std::endl;
  std::cout << pc[1].toString() << std::endl;
  std::cout << pc[pc.size() - 1].toString() << std::endl;

  //create a LJ interaction and set its cutoff
  LennardJonesInteraction lj = LennardJonesInteraction();
  lj.setCutoff(2.5);

  FullPairIterator it = FullPairIterator(pc.size());

  //compute the total potential energy of the system
  compute<LennardJonesInteraction, FullPairIterator>(pc, lj, it, bc);

  return 0;
}
