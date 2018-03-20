/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Tabulated.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "InterpolationLinear.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "python.hpp"

namespace espressopp {
namespace interaction {

void Tabulated::setFilename(int itype, const char* _filename) {
  boost::mpi::communicator world;
  filename = _filename;

  if (itype == 1) {  // create a new InterpolationLinear
    table = make_shared<InterpolationLinear>();
    table->read(world, _filename);
  }

  else if (itype == 2) {  // create a new InterpolationAkima
    table = make_shared<InterpolationAkima>();
    table->read(world, _filename);
  }

  else if (itype == 3) {  // create a new InterpolationCubic
    table = make_shared<InterpolationCubic>();
    table->read(world, _filename);
  }
}

typedef class VerletListInteractionTemplate<Tabulated> VerletListTabulated;
typedef class VerletListAdressInteractionTemplate<Tabulated, Tabulated> VerletListAdressTabulated;
typedef class VerletListHadressInteractionTemplate<Tabulated, Tabulated> VerletListHadressTabulated;
typedef class CellListAllPairsInteractionTemplate<Tabulated> CellListTabulated;
typedef class FixedPairListInteractionTemplate<Tabulated> FixedPairListTabulated;
typedef class FixedPairListTypesInteractionTemplate<Tabulated> FixedPairListTypesTabulated;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void Tabulated::registerPython() {
  using namespace espressopp::python;

  class_<Tabulated, bases<Potential> >("interaction_Tabulated", init<int, const char*, real>())
      .add_property("filename", &Tabulated::getFilename, &Tabulated::setFilename)
      .def_pickle(Tabulated_pickle());

  class_<VerletListTabulated, bases<Interaction> >("interaction_VerletListTabulated",
                                                   init<shared_ptr<VerletList> >())
      .def("setPotential", &VerletListTabulated::setPotential)
      .def("getPotential", &VerletListTabulated::getPotentialPtr);

  class_<VerletListAdressTabulated, bases<Interaction> >(
      "interaction_VerletListAdressTabulated",
      init<shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
      .def("setPotentialAT", &VerletListAdressTabulated::setPotentialAT)
      .def("setPotentialCG", &VerletListAdressTabulated::setPotentialCG);
  ;

  class_<VerletListHadressTabulated, bases<Interaction> >(
      "interaction_VerletListHadressTabulated",
      init<shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
      .def("setPotentialAT", &VerletListHadressTabulated::setPotentialAT)
      .def("setPotentialCG", &VerletListHadressTabulated::setPotentialCG);
  ;

  class_<CellListTabulated, bases<Interaction> >("interaction_CellListTabulated",
                                                 init<shared_ptr<storage::Storage> >())
      .def("setPotential", &CellListTabulated::setPotential);
  ;

  class_<FixedPairListTabulated, bases<Interaction> >(
      "interaction_FixedPairListTabulated",
      init<shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Tabulated> >())
      .def(init<shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Tabulated> >())
      .def("setPotential", &FixedPairListTabulated::setPotential)
      .def("setFixedPairList", &FixedPairListTabulated::setFixedPairList)
      .def("getFixedPairList", &FixedPairListTabulated::getFixedPairList);

  class_<FixedPairListTypesTabulated, bases<Interaction> >(
      "interaction_FixedPairListTypesTabulated",
      init<shared_ptr<System>, shared_ptr<FixedPairList> >())
      .def(init<shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
      .def("setPotential", &FixedPairListTypesTabulated::setPotential)
      .def("getPotential", &FixedPairListTypesTabulated::getPotentialPtr)
      .def("setFixedPairList", &FixedPairListTypesTabulated::setFixedPairList)
      .def("getFixedPairList", &FixedPairListTypesTabulated::getFixedPairList);
}
}
}
