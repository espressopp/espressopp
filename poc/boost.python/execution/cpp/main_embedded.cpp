// the main routine for embedding python
#include <boost/python.hpp>
#include <iostream>

// include the various module python bindings
#include "hello.hpp"

using namespace boost;
using namespace std;

// prototype for the function that python calls when loading a dll
// this function is actually defined by the BOOST_PYTHON_MODULE
// macro. But Python requires exactly this name, so we can be sure
// what the macro does
extern "C" {
  void init_espresso();
}

// call functions that should be called before python starts.
// of course only works in embedded mode
void callBeforePythonHooks(int &argc, char **&argv)
{
  hello::initBeforePython(argc, argv);
}

// list of module(s) that this binary provides
// simply provides the name and the name of the init-routine
// the python routines that use _inittab do not import it
// const, therefore, we cannot put "_espresso" directly
static char espresso_name[] = "_espresso";
static struct _inittab libs[] = {
  {espresso_name, init_espresso},
  {0, 0}
};

// the program entry point that will call python
int main(int argc, char **argv)
{
  callBeforePythonHooks(argc, argv);

  // register the modules that are compiled into this binary
  if (PyImport_ExtendInittab(libs) == -1) {
    cerr << "could not append mpi_test module to the list of loadable modules"
	 << endl;
    exit(-1);
  }

  // and fire up python
  return Py_Main(argc, argv);
}
