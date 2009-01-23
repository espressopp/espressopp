#define RAND_MAX
#define NUMBER_OF_PARTICLES 33000
#define BOX_SIZE 200

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <sys/times.h>

#include "Particle.hpp"
#include "LennardJonesInteraction.hpp"
#include "FullListIterator.hpp"

using namespace std;

#define MAX 4096

real compute (std::vector<Particle> *pc, const Interaction *lj, 
              FullListIterator *it)

{
  real rsq;
  real en;

  en = 0.0;

  int  inum = it->inum;
  int* ilist = it->ilist;

  int*  numneigh   = it->numneigh;
  int** firstneigh = it->firstneigh;

  for (int ii = 0; ii < inum; ii++) {

    int i      = ilist[ii];
    int *jlist = firstneigh[i];
    int jnum   = numneigh[i];

    Particle* Pi = &(*pc)[i];

    real x = Pi->getx();
    real y = Pi->gety();
    real z = Pi->getz();

    for (int jj = 0; jj < jnum; jj++) {

      int j = jlist[jj];

      Particle* Pj = &(*pc)[j];

      real tmp;
      tmp = x - Pj->getx();
      rsq   = tmp * tmp;
      tmp = y - Pj->gety();
      rsq   = tmp * tmp;
      tmp = z - Pj->getz();
      rsq   = tmp * tmp;
      en   += lj->computeEnergy(rsq);
    }

  }

  //write out the total LJ energy

  return en;

}

#define NITER 25

void run(int size) {

  //variables for storing random numbers

  real rx;
  real ry;
  real rz;

  // create a vector to store the particles

  std::vector<Particle> pc(size);

  // assign random positions to the particles on r[0, BOX_SIZE]

  for(int i = 0; i < pc.size(); i++) {
    rx = BOX_SIZE * real(rand()) / RAND_MAX;
    ry = BOX_SIZE * real(rand()) / RAND_MAX;
    rz = BOX_SIZE * real(rand()) / RAND_MAX;
    pc[i] = Particle(rx, ry, rz);
  }

  //print a particle to standard output as a test

  std::cout << pc[1].toString() << std::endl;

  //create a LJ interaction and set its cutoff

  LennardJonesInteraction lj = LennardJonesInteraction();
  lj.setCutoff(2.5);

  // build an iterator for all pairs (i, j) with i < j 

  FullListIterator it = FullListIterator(pc.size());

  real energy;
  int  count = 0;

  for (int iter = 0; iter < NITER; iter++) {

     std::cout << "iter " << iter << " of " << NITER << std::endl;

     energy = compute(&pc, &lj, &it);

     count++;

     // print energy every 5 iterations

     if (count == 5) {
        std::cout << "energy = " << energy << std::endl;
        count = 0;
     }
  }

}

int main() {

  struct tms start_t, stop_t;
  const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

  int sizes[] = {  1000, 1000, 8000 };

  int n = sizeof(sizes) / sizeof(int);

  for (int i = 0; i < n; i++) {

     int size = sizes[i];

     cout << "LJ interaction for " << size << " particles" << endl;

     times(&start_t);
     run (size);
     times(&stop_t);
 
     cout << "LJ interaction for " << size << " particles finished (" << 
              NITER << " iterations)" << endl;
     cout << "user time = "
          << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
          << "s\tsystem time = "
          << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
          << "s" << endl;
     }
}

