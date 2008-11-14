#include "HelloWorld.hpp"

#include <sstream>

#ifdef HAVE_MPI
#include <boost/mpi.hpp>
#endif

#include <boost/python.hpp>

using namespace std;

namespace hello {
  LOG4ESPP_LOGGER(logger, "hello.HelloWorld");

  //////////////////////////////////////////////////
  // IMPLEMENTATION
  //////////////////////////////////////////////////
  void
  HelloWorld::createMessage() {
    string msg = "Hello World!";

#ifdef HAVE_MPI
    boost::mpi::communicator world;
    ostringstream ost;
    ost << "MPI process #" << world.rank() << ": " << msg;

    LOG4ESPP_INFO(logger, pmi::printWorkerId() << "Creating message.");

    // gather messages
    if (pmi::isController()) {
      boost::mpi::gather(world, ost.str(), allMessages, 0);
    } else {
      boost::mpi::gather(world, ost.str(), 0);
    }
#else
    allMessages.push_back(msg);
#endif
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

#ifdef HAVE_MPI
  //////////////////////////////////////////////////
  // REGISTRATION WITH PMI
  //////////////////////////////////////////////////
  // here, you need to register the SERIAL class
  PMI_REGISTER_CLASS(hello::HelloWorld, "hello::HelloWorld");
  PMI_REGISTER_METHOD(hello::HelloWorld, createMessage, "createMessage");
#endif 

#ifdef HAVE_PYTHON  
  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  // here, register the parallel class
  void 
  PHelloWorld::registerPython() {
    using namespace boost::python;
    
#ifdef HAVE_MPI
    new boost::mpi::environment;
#endif

    class_<PHelloWorld>("hello_HelloWorld", init<>())
      .def("getMessages", &PHelloWorld::getMessages)
      .def("__str__", &PHelloWorld::getMessages);
  }
#endif
}
