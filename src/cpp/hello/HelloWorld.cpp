#include "HelloWorld.hpp"
#include "logging/log4espp.hpp"

#include <sstream>

#ifdef HAVE_MPI
#include <boost/mpi.hpp>
#endif

#include <boost/python.hpp>
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
#ifdef HAVE_MPI
      if (pmi::isController())
	pclass.invoke<&HelloWorld::getMessageWorker>();
#endif

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
	BOOST_FOREACH(const string s, allMessages)
	  {	
	    msg += s;
	    msg += "\n";
	  }
#endif

      return msg;
    }


#ifdef HAVE_MPI
    void HelloWorld::getMessageWorker() { getMessage(); }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PMI
    //////////////////////////////////////////////////
    // here, you need to register the SERIAL class
    PMI_REGISTER_CLASS("espresso::hello::HelloWorld", espresso::hello::HelloWorld);
    PMI_REGISTER_METHOD("getMessageWorker", espresso::hello::HelloWorld, getMessageWorker);
#endif 

#ifdef HAVE_PYTHON  
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    // here, register the parallel class
    void 
    HelloWorld::registerPython() {
      using namespace boost::python;
    
      class_<HelloWorld>("hello_HelloWorld", init<>())
	.def("getMessage", &HelloWorld::getMessage)
	.def("__str__", &HelloWorld::getMessage);
    }
#endif
  }
}
