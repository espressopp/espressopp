#include <boost/python.hpp>
#include <vector>
#include <string>
#include <cstring>
#include "mpi_test_internal.hpp"

using namespace std;

namespace python = boost::python;

// the MPI test environment
MPITestInternal *theMPI = 0;

// MPI_Initialize, low level
void initMPI();

// send a command code to all nodes from python
void testMPI(int code);

// MPI_Finalize, called automatically on exit
void finalizeMPI();

void initMPI(int &argc, char **&argv, const char *id)
{
  if (!theMPI)
    theMPI = new MPITestInternal(argc, argv, id);
}

void testMPI(int code)
{
  theMPI->issue_test(code);
}

void finalizeMPI()
{
  // try graceful exit of slaves
  theMPI->issue_test(42);
  delete theMPI; theMPI = 0;
}

/*
 * the entry point for the library version.
 */
BOOST_PYTHON_MODULE(mpi_test)
{
  // similar to boost.python: 
  // convert python argv back to C-style
  python::object sys = python::object(python::handle<>(PyImport_ImportModule("sys")));
  python::list args = python::extract<python::list>(sys.attr("argv"));
  int    argc = python::extract<int>(args.attr("__len__")());
  char **argv = new char *[argc];
  for (int arg = 0; arg < argc; ++arg) {
    argv[arg] = strdup(python::extract<const char*>(args[arg]));
  }

  initMPI(argc, argv, "python module");

  // hand back modified argc and argv to Python
  PySys_SetArgv(argc, argv);

  if (theMPI->slaveMainLoop()) {
    // slave finished, exit - the interpreter must not get control
    // on a slave again
    delete theMPI;
    Py_Exit(1);
  }

  // here, we are the master
  // make sure MPI is properly finalized on exit
  Py_AtExit(finalizeMPI);
  // and add the test command
  python::def("test", testMPI, python::args("code"), "not documented");
}

/*
 * the entry point for the embedded interpreter version
 */
int main(int argc, char **argv)
{
  // MPI initialization outside Python
  initMPI(argc, argv, "python embedded");

  if (theMPI->slaveMainLoop()) {
    // slave finished, exit - we never start Python
    delete theMPI;
    exit(1);
  }

  // here, we are the master

  // list of modules that this binary provides
  // simply provides the name and the name of the init-routine
  char mpi_name[] = "cpp.mpi_test";
  struct _inittab libs[] = {
    {mpi_name, initmpi_test},
    {0, 0}
  };

  // register the modules that are compiled into this binary
  if (PyImport_ExtendInittab(libs) == -1) {
    cerr << "could not append mpi_test module to the list of loadable modules" << endl;
    exit(-1);
  }
  // and fire up normal Python
  return Py_Main(argc, argv);
}
