#include "python.hpp"
#include "System.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "esutil/RNG.hpp"
#include "mpi.hpp"
#include "../contrib/mpi4py/mpi4py-1.2.1/src/include/mpi4py/mpi4py.h"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
# define VT_ON()
# define VT_OFF()
#endif

namespace espresso {

  System::System() {
    comm = mpiWorld;
    CommunicatorIsInitialized = false;
  }

  System::System(python::object _pyobj) {

    //This is defined in mpi4py.h
    //struct __pyx_obj_6mpi4py_3MPI_Comm {
    //  PyObject_HEAD
    //  MPI_Comm ob_mpi;
    //  int flags;
    //};
    

    // Following is some extreme typecasting which we need to convert the
    // pmi python object pmi._MPIcomm into a shared_ptr< boost::mpi::communicator >
    // I have not yet figured out how to do it in a more elegant way
    PyObject *pyobj = _pyobj.ptr();
    // in mpi4py.1.2.1 this has to be:
    __pyx_obj_6mpi4py_3MPI_Comm * pyMPIComm = (__pyx_obj_6mpi4py_3MPI_Comm *) pyobj;
    // in mpi4py.1.2.2 this has to be:
    // PyMPICommObject* pyMPIComm = (PyMPICommObject*) pyobj;
    MPI_Comm * comm_p = &pyMPIComm->ob_mpi;
    shared_ptr< mpi::communicator > newcomm = make_shared< mpi::communicator >(*comm_p, mpi::comm_attach);

    comm = newcomm;
  }

  void System::addInteraction(shared_ptr< interaction::Interaction > ia)
  {
    shortRangeInteractions.push_back(ia);
  }

  void System::removeInteraction(int i)
  {
    shortRangeInteractions.erase(shortRangeInteractions.begin()+i);
  }

  shared_ptr< interaction::Interaction > System::getInteraction(int i)
  {
    return shortRangeInteractions[i];
  }

  int System::getNumberOfInteractions()
  {
	return shortRangeInteractions.size();
  }

  // Scale all coordinates of the system
  void System::scaleVolume(real s) {
	storage->scaleVolume(s);
	bc->scaleVolume(s);
  }

  static void setTrace(bool flag)
  {
     if (flag)
     {
        VT_ON();
     }
     else
     {
        VT_OFF();
     }
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  System::registerPython() {
    using namespace espresso::python;

    class_< System >
      ("System", init<>())
      .def(init< python::object >())
      .def_readwrite("storage", &System::storage)
      .def_readwrite("bc", &System::bc)
      .def_readwrite("rng", &System::rng)
//      .def_readwrite("shortRangeInteractions",
//		     &System::shortRangeInteractions)
      .def_readwrite("skin", &System::skin)
      .def("addInteraction", &System::addInteraction)
      .def("removeInteraction", &System::removeInteraction)
      .def("getInteraction", &System::getInteraction)
      .def("getNumberOfInteractions", &System::getNumberOfInteractions)
      .def("scaleVolume", &System::scaleVolume)
      ;

    def("setTrace", setTrace);
  }
}
