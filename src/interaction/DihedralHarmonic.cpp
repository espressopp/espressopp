/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013,2015 (H)
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

#include "python.hpp"
#include "DihedralHarmonic.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"
#include "FixedQuadrupleListTypesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralHarmonic::registerPython() {
      using namespace espressopp::python;

      class_ <DihedralHarmonic, bases <DihedralPotential> >
      ("interaction_DihedralHarmonic", init< real, real >())
        .add_property("K", &DihedralHarmonic::getK, &DihedralHarmonic::setK)
        .add_property("phi0", &DihedralHarmonic::getPhi0, &DihedralHarmonic::setPhi0)
      ;

      typedef class FixedQuadrupleListInteractionTemplate <DihedralHarmonic>
      FixedQuadrupleListDihedralHarmonic;
      class_ <FixedQuadrupleListDihedralHarmonic, bases <Interaction> >
        ("interaction_FixedQuadrupleListDihedralHarmonic",
                  init< shared_ptr<System>,
                        shared_ptr<FixedQuadrupleList>,
                        shared_ptr<DihedralHarmonic> >())
        .def("setPotential", &FixedQuadrupleListDihedralHarmonic::setPotential)
        .def("getFixedQuadrupleList", &FixedQuadrupleListDihedralHarmonic::getFixedQuadrupleList)
        ;
      typedef class FixedQuadrupleListTypesInteractionTemplate<DihedralHarmonic>
        FixedQuadrupleListTypesDihedralHarmonic;
      class_< FixedQuadrupleListTypesDihedralHarmonic, bases< Interaction > >
        ("interaction_FixedQuadrupleListTypesDihedralHarmonic",
         init< shared_ptr<System>, shared_ptr<FixedQuadrupleList> >())
        .def("setPotential", &FixedQuadrupleListTypesDihedralHarmonic::setPotential)
        .def("getPotential", &FixedQuadrupleListTypesDihedralHarmonic::getPotentialPtr)
        .def("setFixedQuadrupleList", &FixedQuadrupleListTypesDihedralHarmonic::setFixedQuadrupleList)
        .def("getFixedQuadrupleList", &FixedQuadrupleListTypesDihedralHarmonic::getFixedQuadrupleList);
    }
  }
}
