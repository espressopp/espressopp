#include "HelloWorld.hpp"
#include "logging/log4espp.hpp"

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
  const string
  HelloWorld::getMessage() {
    string msg = "Hello World!";

#ifdef HAVE_MPI
    boost::mpi::communicator world;
    ostringstream ost;
    vector<string> allMessages;

    ost << "MPI process #" << world.rank() << ": " << msg;

    LOG4ESPP_INFO(logger, \
		  pmi::printWorkerId() << "Creating message.");

    // gather messages from the tasks
    boost::mpi::gather(world, ost.str(), allMessages, 0);

    msg = "";

    if (pmi::isController())
      // composite message
      for (vector<string>::iterator it = allMessages.begin(); 
	   it != allMessages.end(); it++) {
	msg += *it;
	msg += "\n";
      }
#endif

    return msg;
  }


#ifdef HAVE_MPI
  //////////////////////////////////////////////////
  // REGISTRATION WITH PMI
  //////////////////////////////////////////////////
  // here, you need to register the SERIAL class
  PMI_REGISTER_CLASS("hello::HelloWorld", hello::HelloWorld);
  PMI_REGISTER_METHOD_SPMD("getMessage", hello::HelloWorld, getMessage, const string);
#endif 

#ifdef HAVE_PYTHON  
  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  // here, register the parallel class
  void 
  HelloWorld::registerPython() {
    using namespace boost::python;
    
    class_<PHelloWorld>("hello_HelloWorld", init<>())
      .def("getMessage", &PHelloWorld::getMessage)
      .def("__str__", &PHelloWorld::getMessage);
  }
#endif
}
