#include "python.hpp"
#include "System.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  void System::addInteraction(shared_ptr< interaction::Interaction > ia)
  {
    shortRangeInteractions.push_back(ia);
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
