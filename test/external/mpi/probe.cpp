#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <esutil/Timer.hpp>

const int NUM_TESTS = 1000000;
const int N = 20;

using namespace std;
using namespace espresso::esutil;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  WallTimer timer;

  // Execute the communication test only if there are more than one
  // tasks
  if (MPI::COMM_WORLD.Get_size() > 1) {
    int rank = MPI::COMM_WORLD.Get_rank();

    timer.reset();
    if (rank == 0) {
      double data[N];
      for (int i = 0; i < N; i++)
	data[i] = drand48();
      
      for (int i = 0; i < NUM_TESTS; i++)
	MPI::COMM_WORLD.Send(&data, N, MPI_DOUBLE, 1, 0);
    } else if (rank == 1) {
      double* data = 0;
      MPI::Status status;
      int count;

      for (int i = 0; i < NUM_TESTS; i++) {
	MPI::COMM_WORLD.Probe(0, 0, status);
	count = status.Get_count(MPI_DOUBLE);
	data = static_cast < double* > (realloc(data, sizeof(double)*count));
	MPI::COMM_WORLD.Recv(data, count, MPI_DOUBLE, 0, 0);
      }
    }
    cout << timer << endl;
  }

  if (MPI::COMM_WORLD.Get_size() > 1) {
    timer.reset();
    int rank = MPI::COMM_WORLD.Get_rank();

    if (rank == 0) {
      double data[N];
      for (int i = 0; i < N; i++)
	data[i] = drand48();
      
      for (int i = 0; i < NUM_TESTS; i++)
	MPI::COMM_WORLD.Send(&data, N, MPI_DOUBLE, 1, 0);
    } else if (rank == 1) {
      double data[20];
      MPI::Status status;

      for (int i = 0; i < NUM_TESTS; i++) {
	MPI::COMM_WORLD.Recv(data, N, MPI_DOUBLE, 0, 0);
      }
    }
    cout << timer << endl;
  }

  MPI::Finalize();
  return 0;
}
