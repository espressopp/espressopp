#include "HelloWorld.hpp"

#include <iostream>
#include <boost/mpi.hpp>
#include <boost/python.hpp>

using namespace std;

namespace hello {
  //////////////////////////////////////////////////
  // IMPLEMENTATION
  //////////////////////////////////////////////////
  // gather messages from all processes
  void 
  WHelloWorld::printMessage() {
    boost::mpi::communicator world;
    // TODO
    if (pmi::isController()) return;

    string msg = "Hello World!";

    ostringstream ost;
    ost << "MPI process #" << world.rank() << ": " << msg;

    gather(world, ost.str(), 0);
  }

  string
  CHelloWorld::printMessage() {
    wHello.invoke<&WHelloWorld::printMessage>();

    boost::mpi::communicator world;
    string msg = "Hello World!";

    ostringstream ost;
    ost << "MPI process #" << world.rank() << ": " << msg;

    vector<string> allMessages;
    boost::mpi::gather(world, ost.str(), allMessages, 0);

    string result;
    for (vector<string>::iterator it = allMessages.begin();
	 it != allMessages.end(); it++) {
      result += *it;
      result += "\n";
    }

    return result;
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PMI
  //////////////////////////////////////////////////
  // here, you need to register the SERIAL class
  PMI_REGISTER_CLASS(hello::WHelloWorld, "hello::HelloWorld");
  PMI_REGISTER_METHOD(hello::WHelloWorld, printMessage, "printMessage");
  
  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  // here, register the parallel class
  void 
  CHelloWorld::registerPython() {
    using namespace boost::python;
    
    class_<CHelloWorld>("hello_HelloWorld", init<>())
      .def("printMessage", &CHelloWorld::printMessage)
      .def("__str__", &CHelloWorld::printMessage);
  }
}
