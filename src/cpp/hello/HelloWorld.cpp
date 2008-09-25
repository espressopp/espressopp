#include "HelloWorld.hpp"

#include <sstream>
#include <boost/mpi.hpp>
#include <boost/python.hpp>

using namespace std;

namespace hello {
  LOG4ESPP_LOGGER(logger, "hello.HelloWorld");

  //////////////////////////////////////////////////
  // IMPLEMENTATION
  //////////////////////////////////////////////////
  void
  HelloWorld::createMessage() {
    boost::mpi::communicator world;

    // COMPUTE
    string msg = "Hello World!";
    ostringstream ost;
    ost << "MPI process #" << world.rank() << ": " << msg;

    LOG4ESPP_INFO(logger, pmi::printWorkerId() << "Creating message.");

    // gather messages
    if (pmi::isController()) {
      boost::mpi::gather(world, ost.str(), allMessages, 0);
    } else {
      boost::mpi::gather(world, ost.str(), 0);
    }
  }

  string
  HelloWorld::getMessages() {
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
  PMI_REGISTER_CLASS(hello::HelloWorld, "hello::HelloWorld");
  PMI_REGISTER_METHOD(hello::HelloWorld, createMessage, "createMessage");
  
  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  // here, register the parallel class
  void 
  PHelloWorld::registerPython() {
    using namespace boost::python;
    
    new boost::mpi::environment;

    class_<PHelloWorld>("hello_HelloWorld", init<>())
      .def("getMessages", &PHelloWorld::getMessages)
      .def("__str__", &PHelloWorld::getMessages);
  }
}
