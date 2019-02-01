/*
  Copyright (C) 2017,2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013,2014,2015
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
#include "Tabulated.hpp"
#include "LennardJones.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListAdressCGInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "VerletListHadressCGInteractionTemplate.hpp"
#include "VerletListPIadressInteractionTemplate.hpp"
#include "VerletListPIadressNoDriftInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "FixedPairListPIadressInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {


    void Tabulated::setFilename(int itype, const char* _filename) {
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

    typedef class VerletListInteractionTemplate <Tabulated> VerletListTabulated;
    typedef class VerletListAdressInteractionTemplate <Tabulated, Tabulated> VerletListAdressTabulated;
    typedef class VerletListAdressCGInteractionTemplate <Tabulated> VerletListAdressCGTabulated;
    typedef class VerletListHadressInteractionTemplate <Tabulated, Tabulated> VerletListHadressTabulated;
    typedef class VerletListHadressCGInteractionTemplate <Tabulated> VerletListHadressCGTabulated;
    typedef class VerletListPIadressInteractionTemplate <Tabulated, Tabulated> VerletListPIadressTabulated;
    typedef class VerletListPIadressInteractionTemplate <Tabulated, LennardJones> VerletListPIadressTabulatedLJ;
    typedef class VerletListPIadressNoDriftInteractionTemplate <Tabulated> VerletListPIadressNoDriftTabulated;
    typedef class CellListAllPairsInteractionTemplate <Tabulated> CellListTabulated;
    typedef class FixedPairListInteractionTemplate <Tabulated> FixedPairListTabulated;
    typedef class FixedPairListTypesInteractionTemplate <Tabulated> FixedPairListTypesTabulated;
    typedef class FixedPairListPIadressInteractionTemplate <Tabulated> FixedPairListPIadressTabulated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Tabulated::registerPython() {
      using namespace espressopp::python;

      class_ <Tabulated, bases <Potential> >
        ("interaction_Tabulated", init <int, const char*, real>())
            .add_property("filename", &Tabulated::getFilename, &Tabulated::setFilename)
            .def_pickle(Tabulated_pickle())
        ;

      class_ <VerletListTabulated, bases <Interaction> >
        ("interaction_VerletListTabulated", init <shared_ptr<VerletList> >())
            .def("setPotential", &VerletListTabulated::setPotential)
            .def("getPotential", &VerletListTabulated::getPotentialPtr)
        ;

      class_< VerletListAdressCGTabulated, bases< Interaction > >
        ("interaction_VerletListAdressCGTabulated",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressCGTabulated::getVerletList)
        .def("setPotential", &VerletListAdressCGTabulated::setPotential)
        .def("getPotential", &VerletListAdressCGTabulated::getPotentialPtr)
      ;

      class_ <VerletListAdressTabulated, bases <Interaction> >
        ("interaction_VerletListAdressTabulated",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress> >()
                )
            .def("setPotentialAT", &VerletListAdressTabulated::setPotentialAT)
            .def("setPotentialCG", &VerletListAdressTabulated::setPotentialCG);
        ;

      class_< VerletListHadressCGTabulated, bases< Interaction > >
        ("interaction_VerletListHadressCGTabulated",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressCGTabulated::getVerletList)
        .def("setPotential", &VerletListHadressCGTabulated::setPotential)
        .def("getPotential", &VerletListHadressCGTabulated::getPotentialPtr)
      ;

      class_ <VerletListHadressTabulated, bases <Interaction> >
        ("interaction_VerletListHadressTabulated",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress> >()
                )
            .def("setPotentialAT", &VerletListHadressTabulated::setPotentialAT)
            .def("setPotentialCG", &VerletListHadressTabulated::setPotentialCG);
        ;

      class_ <VerletListPIadressTabulated, bases <Interaction> >
        ("interaction_VerletListPIadressTabulated",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress>,
                 int,
                 bool>()
                )
            .def("setPotentialQM", &VerletListPIadressTabulated::setPotentialQM)
            .def("setPotentialCL", &VerletListPIadressTabulated::setPotentialCL)
            .def("setVerletList", &VerletListPIadressTabulated::setVerletList)
            .def("getVerletList", &VerletListPIadressTabulated::getVerletList)
            .def("setFixedTupleList", &VerletListPIadressTabulated::setFixedTupleList)
            .def("getFixedTupleList", &VerletListPIadressTabulated::getFixedTupleList)
            .def("setNTrotter", &VerletListPIadressTabulated::setNTrotter)
            .def("getNTrotter", &VerletListPIadressTabulated::getNTrotter)
            .def("setSpeedup", &VerletListPIadressTabulated::setSpeedup)
            .def("getSpeedup", &VerletListPIadressTabulated::getSpeedup);
        ;

      class_ <VerletListPIadressTabulatedLJ, bases <Interaction> >
        ("interaction_VerletListPIadressTabulatedLJ",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress>,
                 int,
                 bool>()
                )
            .def("setPotentialQM", &VerletListPIadressTabulatedLJ::setPotentialQM)
            .def("setPotentialCL", &VerletListPIadressTabulatedLJ::setPotentialCL)
            .def("setVerletList", &VerletListPIadressTabulatedLJ::setVerletList)
            .def("getVerletList", &VerletListPIadressTabulatedLJ::getVerletList)
            .def("setFixedTupleList", &VerletListPIadressTabulatedLJ::setFixedTupleList)
            .def("getFixedTupleList", &VerletListPIadressTabulatedLJ::getFixedTupleList)
            .def("setNTrotter", &VerletListPIadressTabulatedLJ::setNTrotter)
            .def("getNTrotter", &VerletListPIadressTabulatedLJ::getNTrotter)
            .def("setSpeedup", &VerletListPIadressTabulatedLJ::setSpeedup)
            .def("getSpeedup", &VerletListPIadressTabulatedLJ::getSpeedup);
        ;

      class_ <VerletListPIadressNoDriftTabulated, bases <Interaction> >
        ("interaction_VerletListPIadressNoDriftTabulated",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress>,
                 int,
                 bool>()
                )
            .def("setPotential", &VerletListPIadressNoDriftTabulated::setPotential)
            .def("setVerletList", &VerletListPIadressNoDriftTabulated::setVerletList)
            .def("getVerletList", &VerletListPIadressNoDriftTabulated::getVerletList)
            .def("setFixedTupleList", &VerletListPIadressNoDriftTabulated::setFixedTupleList)
            .def("getFixedTupleList", &VerletListPIadressNoDriftTabulated::getFixedTupleList)
            .def("setNTrotter", &VerletListPIadressNoDriftTabulated::setNTrotter)
            .def("getNTrotter", &VerletListPIadressNoDriftTabulated::getNTrotter)
            .def("setSpeedup", &VerletListPIadressNoDriftTabulated::setSpeedup)
            .def("getSpeedup", &VerletListPIadressNoDriftTabulated::getSpeedup);
        ;

      class_ <CellListTabulated, bases <Interaction> >
        ("interaction_CellListTabulated", init <shared_ptr <storage::Storage> >())
            .def("setPotential", &CellListTabulated::setPotential);
        ;

      class_ <FixedPairListTabulated, bases <Interaction> >
        ("interaction_FixedPairListTabulated",
          init <shared_ptr<System>,
                shared_ptr<FixedPairList>,
                shared_ptr<Tabulated> >()
        )
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Tabulated> >())
        .def("setPotential", &FixedPairListTabulated::setPotential)
        .def("setFixedPairList", &FixedPairListTabulated::setFixedPairList)
        .def("getFixedPairList", &FixedPairListTabulated::getFixedPairList);

      class_< FixedPairListTypesTabulated, bases< Interaction > >
          ("interaction_FixedPairListTypesTabulated",
           init< shared_ptr<System>, shared_ptr<FixedPairList> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
          .def("setPotential", &FixedPairListTypesTabulated::setPotential)
          .def("getPotential", &FixedPairListTypesTabulated::getPotentialPtr)
          .def("setFixedPairList", &FixedPairListTypesTabulated::setFixedPairList)
          .def("getFixedPairList", &FixedPairListTypesTabulated::getFixedPairList);

      class_ <FixedPairListPIadressTabulated, bases <Interaction> >
        ("interaction_FixedPairListPIadressTabulated",
          init <shared_ptr<System>,
                shared_ptr<FixedPairList>,
                shared_ptr<FixedTupleListAdress>,
                shared_ptr<Tabulated>,
                int,
                bool>()
        )
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<FixedTupleListAdress>, shared_ptr<Tabulated>, int, bool>())
        .def("setPotential", &FixedPairListPIadressTabulated::setPotential)
        .def("getPotential", &FixedPairListPIadressTabulated::getPotential)
        .def("setFixedPairList", &FixedPairListPIadressTabulated::setFixedPairList)
        .def("getFixedPairList", &FixedPairListPIadressTabulated::getFixedPairList)
        .def("setFixedTupleList", &FixedPairListPIadressTabulated::setFixedTupleList)
        .def("getFixedTupleList", &FixedPairListPIadressTabulated::getFixedTupleList)
        .def("setNTrotter", &FixedPairListPIadressTabulated::setNTrotter)
        .def("getNTrotter", &FixedPairListPIadressTabulated::getNTrotter)
        .def("setSpeedup", &FixedPairListPIadressTabulated::setSpeedup)
        .def("getSpeedup", &FixedPairListPIadressTabulated::getSpeedup);
        ;
    }

  }
}
