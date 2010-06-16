#include "python.hpp"
#include "System.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "esutil/RNG.hpp"
#include "mpi.hpp"

namespace espresso {

  void System::addInteraction(shared_ptr< interaction::Interaction > ia)
  {
    shortRangeInteractions.push_back(ia);
  }

  int System::sum(int sumLocal)
  {
    int sumGlobal;
    boost::mpi::reduce(*comm, sumLocal, sumGlobal, std::plus<int>(), 0);
    return sumGlobal;
  }

  real System::sum(real sumLocal)
  {
    real sumGlobal;
    boost::mpi::reduce(*comm, sumLocal, sumGlobal, std::plus<real>(), 0);
    return sumGlobal;
  }

  void System::sum(real* sumLocal, real* sumGlobal, int N)
  {
    boost::mpi::reduce(*comm, sumLocal, N, sumGlobal, std::plus<real>(), 0);
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  System::registerPython() {
    using namespace espresso::python;

    class_< System >("System")
      .def_readwrite("storage", &System::storage)
      .def_readwrite("bc", &System::bc)
      .def_readwrite("rng", &System::rng)
      .def_readwrite("shortRangeInteractions", 
		     &System::shortRangeInteractions)
      .def_readwrite("skin", &System::skin)
      .def_readwrite("comm", &System::comm)
      .def("addInteraction", &System::addInteraction)
      ;
  }
}
