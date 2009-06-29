#include <boost/mpi.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

using namespace std;
namespace mpi = boost::mpi;

class Hello {
public:
  void printMessage(ostream &out = cout) {
    out << "Hello World!";
  }
  void doit(ostream &out = cout) {
    out << "Do it!";
  }
};

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  Hello hello;
  void (Hello::*method_ptr)(ostream &out) = &Hello::doit;

  stringstream msg;
  msg << "Worker " << world.rank() << ": ";
  (hello.*method_ptr)(msg);

  // how to get a print of the pointer
  union {
    void (Hello::*method_ptr)(ostream &out);
    unsigned char *cp;
  } u = { method_ptr };

  msg << " (";
  msg.fill('0');
  msg << hex;
  for (unsigned idx = sizeof(method_ptr); idx > 0; idx--)
    msg << setw(2) << static_cast<unsigned>(u.cp[idx]);
  msg << ")";

  if (world.rank() == 0) {
    vector<string> all_messages;
    gather(world, msg.str(), all_messages, 0);
    for (int proc = 0; proc < world.size(); ++proc)
      cout << all_messages[proc] << endl;
  } else {
    gather(world, msg.str(), 0);
  }

  return 0;
} 

