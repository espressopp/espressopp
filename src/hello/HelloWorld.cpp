#include "HelloWorld.hpp"

#include <sstream>

#include <logging.hpp>
#include <mpi.hpp>
#include <python.hpp>
#include <boost/foreach.hpp>

using namespace std;

namespace espresso {
  namespace hello {
    LOG4ESPP_LOGGER(logger, "hello.HelloWorld");
    
    //////////////////////////////////////////////////
    // IMPLEMENTATION
    //////////////////////////////////////////////////
    const string
    HelloWorld::getMessage() {
      // collect the messages of the processors
      ostringstream ost;
      boost::mpi::communicator world;
      ost << "MPI process #" << world.rank() << ": Hello World!";

      LOG4ESPP_INFO(logger,			\
		    "Created message: "		\
		    << ost.str());

      return ost.str();
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    HelloWorld::registerPython() {
      using namespace boost::python;
    
      class_<HelloWorld>("hello_HelloWorld", init<>())
	.def("getMessage", &HelloWorld::getMessage)
	.def("__str__", &HelloWorld::getMessage);
    }
  }
}
