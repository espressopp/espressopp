// Example for using PMI and a proxy class.
//
// Demonstrates, how a parallel proxy class can be created for a
// class, that can be called as though it was the original object.

#include <mpi.h>
#include <iostream>
#include <pmi/pmi.hpp>

using namespace std;

// define a very simple class
class HelloWorld {
public:
  void printMessage() {
    cout << "Worker " << pmi::getWorkerId() 
	 << ": Hello World!" 
	 << endl;
  }
};

// register the class and method with PMI
PMI_REGISTER_CLASS(HelloWorld, "HelloWorld");
PMI_REGISTER_METHOD(HelloWorld, printMessage, "printMessage");

// define the parallel proxy class
class PHelloWorld 
  : public pmi::ParallelObject<HelloWorld> {
public:
  // define its methods
  PMI_PROXY_METHOD(printMessage);
};

int main(int argc, char* argv[]) {
  // Required by MPI
  MPI::Init(argc, argv);
  // Required by the logging system
  LOG4ESPP_CONFIGURE();

  if (!pmi::mainLoop()) {
    // Create an instance of the proxy class
    PHelloWorld hello;
    
    // Call the method
    hello.printMessage();

    // Stop the workers
    endWorkers();
  }
    
  return 0;
}

