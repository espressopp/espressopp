#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  MPI::Init(argc, argv);

  if (MPI::COMM_WORLD.Get_size() > 1) {
    int rank = MPI::COMM_WORLD.Get_rank();

    int data;

    // first test: send data before barrier, receive after barrier
    if (rank == 0) {
      data = 42;
      MPI::COMM_WORLD.Isend(&data, 1, MPI_INT, 1, 0);
      cout << "Task 0: Sent " << data << " before barrier!" << endl;
      MPI::COMM_WORLD.Barrier();
    } else if (rank == 1) {
      MPI::COMM_WORLD.Barrier();
      MPI::Request request = MPI::COMM_WORLD.Irecv(&data, 1, MPI_INT, 0, 0);
      if (request.Test()) {
	cout << "Task 1: Received " << data << " after barrier!" << endl;
      } else {
	cout << "Task 1: Did not receive data after barrier!" << endl;
      }
    }

    // second test: don't send anything
    if (rank == 0) {
      MPI::COMM_WORLD.Barrier();
      cout << "Task 0: Did not send and data before barrier!" << endl;
    } else if (rank == 1) {
      MPI::COMM_WORLD.Barrier();
      MPI::Request request = MPI::COMM_WORLD.Irecv(&data, 1, MPI_INT, 0, 0);
      if (request.Test()) {
	cout << "Task 1: Received " << data << " after barrier!" << endl;
      } else {
	cout << "Task 1: Did not receive data after barrier!" << endl;
      }
    }

    MPI::Finalize();
    return 0;
  } else {
    cout << "Please start at least two MPI tasks!" << endl;
    MPI::Finalize();
    return 1;
  }
}
