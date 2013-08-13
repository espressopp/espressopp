#include "python.hpp"
#include "Tabulated.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
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
    typedef class VerletListHadressInteractionTemplate <Tabulated, Tabulated> VerletListHadressTabulated;
    typedef class CellListAllPairsInteractionTemplate <Tabulated> CellListTabulated;
    typedef class FixedPairListInteractionTemplate <Tabulated> FixedPairListTabulated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Tabulated::registerPython() {
      using namespace espresso::python;
     
      class_ <Tabulated, bases <Potential> >
        ("interaction_Tabulated", init <int, const char*, real>())
            .add_property("filename", &Tabulated::getFilename, &Tabulated::setFilename)
            .def_pickle(Tabulated_pickle())
        ;
     
      class_ <VerletListTabulated, bases <Interaction> > 
        ("interaction_VerletListTabulated", init <shared_ptr<VerletList> >())
            .def("setPotential", &VerletListTabulated::setPotential, return_value_policy< reference_existing_object >())
            .def("getPotential", &VerletListTabulated::getPotential, return_value_policy< reference_existing_object >())
            .def("clonePotential", &VerletListTabulated::clonePotential)
        ;

      class_ <VerletListAdressTabulated, bases <Interaction> >
        ("interaction_VerletListAdressTabulated",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress> >()
                )
            .def("setPotentialAT", &VerletListAdressTabulated::setPotentialAT)
            .def("setPotentialCG", &VerletListAdressTabulated::setPotentialCG);
        ;
        
      class_ <VerletListHadressTabulated, bases <Interaction> >
        ("interaction_VerletListHadressTabulated",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress> >()
                )
            .def("setPotentialAT", &VerletListHadressTabulated::setPotentialAT)
            .def("setPotentialCG", &VerletListHadressTabulated::setPotentialCG);
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
        ;
        
        
    }
    
  }
}
