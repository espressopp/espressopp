#include "HelloWorld.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <boost/mpi.hpp>
#include <boost/python.hpp>

using namespace std;

namespace hello {
  //////////////////////////////////////////////////
  // IMPLEMENTATION
  //////////////////////////////////////////////////
  // gather messages from all processes
  void 
  HelloWorld::print() {
    boost::mpi::communicator world;
    string msg = "Hello World!";

    ostringstream ost;
    ost << "MPI process #" << world.rank() << ": " << msg;

    if (world.rank() == 0) {
      // gather the messages from all processors
      vector<string> allMessages;
      boost::mpi::gather(world, ost.str(), allMessages, 0);

      // write out all messages
      for (vector<string>::iterator it = allMessages.begin();
	   it != allMessages.end(); it++)
	cout << *it << endl;
    } else {
      gather(world, ost.str(), 0);
    }
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PMI
  //////////////////////////////////////////////////
  // here, you need to register the SERIAL class
  PMI_REGISTER_CLASS(hello::HelloWorld, "hello::HelloWorld");
  PMI_REGISTER_METHOD(hello::HelloWorld, print, "print");

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  // here, register the parallel class
  void 
  PHelloWorld::registerPython() {
    using namespace boost::python;
    
    class_<PHelloWorld>("hello_HelloWorld", init<>())
      .def("print", &PHelloWorld::print)
      .def("__str__", &PHelloWorld::print);
  }
}
