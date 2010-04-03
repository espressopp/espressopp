#include <main/espresso_common.hpp>
#include <mpi.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace boost;
using namespace espresso;

void parent()
{
  mpi::packed_oarchive send(*mpiWorld);
  int out;
  out = 1;
  send << out;
  out = 2;
  send << out;
  out = 3;
  send << out;

  mpi::broadcast(*mpiWorld, send, 0);
}

void son()
{
  mpi::packed_iarchive recv(*mpiWorld);
  int tmp;

  mpi::broadcast(*mpiWorld, recv, 0);
  
  recv >> tmp;
  cout << mpiWorld->rank() << ": 1 = " << tmp << endl;

  recv >> tmp;
  cout << mpiWorld->rank() << ": 2 = " << tmp << endl;

  recv >> tmp;
  cout << mpiWorld->rank() << ": 3 = " << tmp << endl;
}

int main()
{
  initMPIEnv();

  if (mpiWorld->rank() == 0) {
    parent();
  }
  else {
    son();
  }
  
  finalizeMPIEnv();
  return 0;	
}
