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
      pmiObject.invoke<&HelloWorld::getMessageWorker>();
      
      string msg = "Hello World!";

      // collect the messages of the processors
      boost::mpi::communicator world;
      ostringstream ost;
      vector<string> allMessages;

      ost << "MPI process #" << world.rank() << ": " << msg;

      LOG4ESPP_INFO(logger, \
		    pmi::printWorkerId() << "Created message: " \
		    << ost.str());

      // gather messages from the tasks
      boost::mpi::gather(world, ost.str(), allMessages, 0);

      msg = "";

      // compose message
      BOOST_FOREACH(const string s, allMessages)
	{	
	  msg += s;
	  msg += "\n";
	}

      return msg;
    }


    void HelloWorld::getMessageWorker() {
      string msg = "Hello World!";

      // send the messages to the controller
      boost::mpi::communicator world;
      ostringstream ost;

      ost << "MPI process #" << world.rank() << ": " << msg;

      LOG4ESPP_INFO(logger, \
		    pmi::printWorkerId() << "Created message: " \
		    << ost.str());

      // gather messages from the tasks
      boost::mpi::gather(world, ost.str(), 0);
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PMI
    //////////////////////////////////////////////////
    PMI_REGISTER_CLASS(espresso::hello::HelloWorld);
    PMI_REGISTER_METHOD(espresso::hello::HelloWorld, getMessageWorker);

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
