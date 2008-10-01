#include <boost/mpi.hpp>
#include <iostream>
#include <boost/serialization/string.hpp>

namespace mpi = boost::mpi;

int main(int argc, char* argv[]) 

{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  int TAG_0 = 0;
  int TAG_1 = 1;

  if (world.size() == 1) {
    std::string test("Hello world");
    broadcast(world, test, 0);
    std::cout << test << "!" << std::endl;
    return 0;
  }

  if (world.rank() == 0) {

    int myPartnerRank = 1;

    world.send(myPartnerRank, TAG_0, std::string("Hello"));
    std::string msg;
    world.recv(myPartnerRank, TAG_1, msg);
    std::cout << msg << "!" << std::endl;
  } else {
 
    int myPartnerRank = 0;

    std::string msg;
    world.recv(myPartnerRank, TAG_0, msg);
    std::cout << msg << ", ";
    std::cout.flush();
    world.send(myPartnerRank, TAG_1, std::string("world"));
  }

  return 0;
}

