#include "mpi_test_internal.hpp"
#include <iostream>

using namespace std;

mpi::environment *MPITestInternal::env = 0;

MPITestInternal::MPITestInternal(int &argc, char **&argv, const char *id)
{
  if (env != 0) {
    cerr << "MPITestInternal instantiated twice" << endl;
    exit(-1);
  }
  env = new mpi::environment(argc, argv);

  mpi::communicator world;
  cout << "rank " << world.rank() << " started from " << id << endl;
}

MPITestInternal::~MPITestInternal()
{
  mpi::communicator world;
  cout << "rank " << world.rank() << " finalizing" << endl;
  delete env; env = 0;
}

bool MPITestInternal::slaveMainLoop()
{
  int code;
  mpi::communicator world;

  if (world.rank() == 0) {
    return false;
  }

  for (;;) {
    broadcast(world, code, 0);
    cout << "rank " << world.rank() << " got job " << code << endl;
    if (code == 42) {
      break;
    }
  }

  return true;
}

void MPITestInternal::issue_test(int code)
{
  mpi::communicator world;

  broadcast(world, code, 0);
  cout << "master sent job " << code << endl;
}
