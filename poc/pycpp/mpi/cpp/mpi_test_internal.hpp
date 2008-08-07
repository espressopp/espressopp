#ifndef INT_HPP
#define INT_HPP

#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

class MPITestInternal {
  class environment;

private:
  static mpi::environment *env;

public:
  MPITestInternal(int &argc, char **&argv, const char *id);
  ~MPITestInternal();

  bool slaveMainLoop();

  void issue_test(int test);
};

#endif
