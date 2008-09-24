// Basic example for using PMI.
//
// Demonstrates the most basic usage of PMI.
//
// Execute via:
//
//   > mpiexec -n 2 hello

#include <iostream>
#include <pmi/pmi.hpp>
#include <pmi/log4espp.hpp>
#include <mpi.h>

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

int main(int argc, char* argv[]) {
  // Required by MPI
  MPI::Init(argc, argv);
  // Required by the logging system
  LOG4ESPP_CONFIGURE();

  // mainLoop will return "false" only on the controller
  if (!pmi::mainLoop()) {
    // Create a parallel instance of the object
    pmi::ParallelObject<HelloWorld> hello;
    
    // Call a method of the object
    hello.invoke<&HelloWorld::printMessage>();

    // Stop the workers
    endWorkers();
  } 

    
  return 0;
}

