/*
  Copyright (C) 2014,2016
      Jakub Krajniak
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

#include "python.hpp"
#include "DihedralHarmonicNCos.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"
#include "FixedQuadrupleListTypesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralHarmonicNCos::registerPython() {
      using namespace espressopp::python;

      class_ <DihedralHarmonicNCos, bases <DihedralPotential> >
      ("interaction_DihedralHarmonicNCos", init< real, real, int >())
        .add_property("K", &DihedralHarmonicNCos::getK, &DihedralHarmonicNCos::setK)
        .add_property("phi0", &DihedralHarmonicNCos::getPhi0, &DihedralHarmonicNCos::setPhi0)
        .add_property("multiplicity", &DihedralHarmonicNCos::getMultiplicity, &DihedralHarmonicNCos::setMultiplicity)
      ;

      typedef class FixedQuadrupleListInteractionTemplate <DihedralHarmonicNCos>
      FixedQuadrupleListDihedralHarmonicNCos;
      class_ <FixedQuadrupleListDihedralHarmonicNCos, bases <Interaction> >
        ("interaction_FixedQuadrupleListDihedralHarmonicNCos",
                  init< shared_ptr<System>,
                        shared_ptr<FixedQuadrupleList>,
                        shared_ptr<DihedralHarmonicNCos> >())
        .def(init< shared_ptr<System>,
                   shared_ptr<FixedQuadrupleListAdress>,
                   shared_ptr<DihedralHarmonicNCos> >())
        .def("setPotential", &FixedQuadrupleListDihedralHarmonicNCos::setPotential)
        .def("getFixedQuadrupleList", &FixedQuadrupleListDihedralHarmonicNCos::getFixedQuadrupleList)
        ;

      typedef class FixedQuadrupleListTypesInteractionTemplate<DihedralHarmonicNCos>
          FixedQuadrupleListTypesDihedralHarmonicNCos;
      class_< FixedQuadrupleListTypesDihedralHarmonicNCos, bases< Interaction > >
          ("interaction_FixedQuadrupleListTypesDihedralHarmonicNCos",
           init< shared_ptr<System>, shared_ptr<FixedQuadrupleList> >())
          .def("setPotential", &FixedQuadrupleListTypesDihedralHarmonicNCos::setPotential)
          .def("getPotential", &FixedQuadrupleListTypesDihedralHarmonicNCos::getPotentialPtr)
          .def("setFixedQuadrupleList", &FixedQuadrupleListTypesDihedralHarmonicNCos::setFixedQuadrupleList)
          .def("getFixedQuadrupleList", &FixedQuadrupleListTypesDihedralHarmonicNCos::getFixedQuadrupleList);
    }
  }
}
