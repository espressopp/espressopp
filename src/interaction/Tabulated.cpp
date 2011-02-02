#include "python.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {


    void Tabulated::setFilename(const char* _filename) {
      boost::mpi::communicator world;
     
      filename = _filename;
     
      // create a new InterpolationTable
      table = make_shared <InterpolationTable> ();
      table->read(world, _filename);
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
        ("interaction_Tabulated", init <const char*, real>())
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
