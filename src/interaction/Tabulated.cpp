#include "python.hpp"
#include "Tabulated.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"
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
    typedef class CellListAllPairsInteractionTemplate <Tabulated> CellListTabulated;
    typedef class FixedPairListInteractionTemplate <Tabulated> FixedPairListTabulated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Tabulated::registerPython() {
      using namespace espresso::python;
     
      class_ <Tabulated, bases <Potential> >
        ("interaction_Tabulated", init <int, const char*, real>())
            .add_property("filename", &Tabulated::getFilename, &Tabulated::setFilename);
     
      class_ <VerletListTabulated, bases <Interaction> > 
        ("interaction_VerletListTabulated", init <shared_ptr<VerletList> >())
            .def("setPotential", &VerletListTabulated::setPotential);
        ;
     
      class_ <CellListTabulated, bases <Interaction> > 
        ("interaction_CellListTabulated", init <shared_ptr <storage::Storage> >())
            .def("setPotential", &CellListTabulated::setPotential);
        ;
        
      class_ <FixedPairListTabulated, bases <Interaction> > 
        ("interaction_FixedPairListTabulated",
            init <shared_ptr<System>, shared_ptr <FixedPairList> >())
            .def("setPotential", &FixedPairListTabulated::setPotential);
        ;
        
        
    }
    
  }
}
