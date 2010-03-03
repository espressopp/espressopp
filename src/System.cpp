#include "System.hpp"
#include <python.hpp>
#include "BC.hpp"
#include "Storage.hpp"
#include "Interaction.hpp"
#include "esutil/RNG.hpp"

namespace espresso {


  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  System::registerPython() {
    using namespace espresso::python;

    class_< System >("System")
      .def_readwrite("name",&System::name)
      .def_readwrite("storage",&System::storage)
      .def_readwrite("bc",&System::bc)
      .def_readwrite("shortRangeInteractions",&System::shortRangeInteractions)
      .def_readwrite("rng",&System::rng)
      .def_readwrite("skin",&System::skin)
      ;
  }
}
