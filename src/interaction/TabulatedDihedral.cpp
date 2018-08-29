/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research
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

#include "python.hpp"
#include "TabulatedDihedral.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"
#include "FixedQuadrupleListTypesInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {

        void TabulatedDihedral::setFilename(int itype, const char* _filename) {
            boost::mpi::communicator world;
            filename = _filename;


            if (itype == 1) { // create a new InterpolationLinear
                table = make_shared <InterpolationLinear> ();
                table->read(world, _filename);
            }

            else if (itype == 2) { // create a new InterpolationAkima
                table = make_shared <InterpolationAkima> ();
                table->read(world, _filename);
            }

            else if (itype == 3) { // create a new InterpolationCubic
                table = make_shared <InterpolationCubic> ();
                table->read(world, _filename);
            }
        }

        typedef class FixedQuadrupleListInteractionTemplate <TabulatedDihedral>
            FixedQuadrupleListTabulatedDihedral;

        typedef class FixedQuadrupleListTypesInteractionTemplate<TabulatedDihedral>
            FixedQuadrupleListTypesTabulatedDihedral;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedDihedral::registerPython() {
            using namespace espressopp::python;

            class_ <TabulatedDihedral, bases <DihedralPotential> >
                ("interaction_TabulatedDihedral", init <int, const char*>())
                .add_property("filename", &TabulatedDihedral::getFilename, &TabulatedDihedral::setFilename)
                .def_pickle(TabulatedDihedral_pickle())
                ;

            class_ <FixedQuadrupleListTabulatedDihedral, bases <Interaction> >
                ("interaction_FixedQuadrupleListTabulatedDihedral",
                        init <shared_ptr<System>,
                              shared_ptr<FixedQuadrupleList>,
                              shared_ptr<TabulatedDihedral> >())
                .def("setPotential", &FixedQuadrupleListTabulatedDihedral::setPotential)
                .def("getFixedQuadrupleList", &FixedQuadrupleListTabulatedDihedral::getFixedQuadrupleList);

            class_< FixedQuadrupleListTypesTabulatedDihedral, bases< Interaction > >
                ("interaction_FixedQuadrupleListTypesTabulatedDihedral",
                 init< shared_ptr<System>, shared_ptr<FixedQuadrupleList> >())
                .def("setPotential", &FixedQuadrupleListTypesTabulatedDihedral::setPotential)
                .def("getPotential", &FixedQuadrupleListTypesTabulatedDihedral::getPotentialPtr)
                .def("setFixedQuadrupleList", &FixedQuadrupleListTypesTabulatedDihedral::setFixedQuadrupleList)
                .def("getFixedQuadrupleList", &FixedQuadrupleListTypesTabulatedDihedral::getFixedQuadrupleList);

        }

    } // ns interaction
} // ns espressopp
