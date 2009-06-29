#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <iostream>

using namespace std;
namespace mpi = boost::mpi;

template < class T >
void
testsend(T msg) {
  mpi::communicator world;
  if (world.rank() == 0) {
    cout << "Sending: " << msg << endl;
    world.send(1, 0, msg);
  } else {
    T received;
    world.recv(0, 0, received);
    cout << "Received: " << received << endl;
  }
}

// Test sending enums
typedef enum {
  A, B, C
} EnumType;

ostream& operator<<(ostream& out, const EnumType &anEnum) {
  out << "EnumType=";
  switch (anEnum) {
  case A: out << "A"; break;
  case B: out << "B"; break;
  case C: out << "C"; break;
  }
  return out;
}

void
send_enum() {
  EnumType msg = A;
  testsend(msg);
}

// Test sending a simple class
class SimpleClass {
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & anEnum;
    ar & anInt;
    ar & anotherInt;
  }

public:
  EnumType anEnum;
  int anInt;
  int anotherInt;
};

ostream& operator<<(ostream& out, const SimpleClass &a) {
  out << "SimpleClass= { " << a.anEnum << ","
      << a.anInt << ","
      << a.anotherInt << " }";
  return out;
}

void
send_SimpleClass() {
  SimpleClass msg;
  msg.anEnum = A;
  msg.anInt = 42;
  msg.anotherInt = 52;
  testsend(msg);
}

int main(int argc, char* argv[]) {
  mpi::environment env(argc, argv);

  send_enum();
  send_SimpleClass();

  return 0;
}
