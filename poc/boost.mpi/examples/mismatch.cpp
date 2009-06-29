#include <boost/mpi.hpp>
#include <iostream>

using namespace std;
namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    int i;
    world.recv(1, 0, i);
    cout << "Received: " << i << endl;

    world.recv(1, 1, i);
    cout << "Received: " << i << endl;
  } else {
    string msg;
    msg = "Hello";
    world.send(0, 1, msg);
    cout << "Sent: " << msg << endl;

    msg = "World!";
    world.send(0, 0, msg);
    cout << "Sent: " << msg << endl;
  }

  return 0;

}
