//  (C) Copyright Thomas Brandes, SCAI Fraunhofer

// MpiTest

#include <mpi.h>
#include <iostream>

using namespace std;
using namespace MPI;

bool testMPIBroadcasting() {
  const int broadcastValue = 13;
  int val = 0;
  int rank = COMM_WORLD.Get_rank();

  if (rank == 0) {
    val = broadcastValue;
  }

  COMM_WORLD.Bcast(&val, 1, MPI_INT, 0);

  // verify that now all processors have the value
  if (val != broadcastValue) {
    cerr << "rank " << rank << ": broadcasting failed (expected to get "
	 << broadcastValue << ", got " << val << ")" << endl;
    return false;
  }

  return true;
}

bool testMPISum() {
  int rank = COMM_WORLD.Get_rank();
  // val is send buffer, sum is receive buffer 
  int val = COMM_WORLD.Get_rank();
  int sum  = 1711;  // must be overwritten

  COMM_WORLD.Allreduce(&val, &sum, 1, MPI_INT, MPI_SUM);

  int size = COMM_WORLD.Get_size();
  // expected output (sum of processor numbers)
  size = size * (size - 1) / 2;

  if (sum  != size) {
    cerr << "rank " << rank << ": all2all sum failed (expected to get "
	 << size << ", got " << sum << ")" << endl;
    return false;
  }
}

int main(int argc, char **argv) {
  Init(argc, argv);

  if (!testMPIBroadcasting() ||
      !testMPISum()) {
    return EXIT_FAILURE;
  }

  Finalize();

  return EXIT_SUCCESS;
}
