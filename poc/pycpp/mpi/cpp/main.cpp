#include "mpi_test_internal.hpp"

int main(int argc, char **argv)
{
  MPITestInternal test(argc, argv, "C++");

  if (!test.slaveMainLoop()) {
    test.issue_test(2);
    test.issue_test(42);
  }
}
