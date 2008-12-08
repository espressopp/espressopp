#include <mpi.h>
#include <iostream>

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  // Execute the communication test only if there are more than one
  // tasks
  if (MPI::COMM_WORLD.Get_size() > 1) {
    int rank = MPI::COMM_WORLD.Get_rank();
    if (rank == 0) {
      int value = 17;
      MPI::COMM_WORLD.Send(&value, 1, MPI_INT, 1, 0);
      std::cout << "Rank 0 OK!" << std::endl;
    } else if (rank == 1) {
      int value;
      MPI::COMM_WORLD.Recv(&value, 1, MPI_INT, 0, 0);
      if (value == 17) std::cout << "Rank 1 OK!" << std::endl;
    } 
  }
  
  MPI::Finalize();
  return 0;
} 
