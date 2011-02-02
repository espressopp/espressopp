#include "python.hpp"
#include "TabulatedDihedral.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"

namespace espresso {
    namespace interaction {
        
        void TabulatedDihedral::setFilename(const char* _filename) {
            boost::mpi::communicator world;
            filename = _filename;
            // create a new InterpolationTable
            table = make_shared <InterpolationTable> ();
            table->read(world, _filename);
        }

        typedef class FixedQuadrupleListInteractionTemplate <TabulatedDihedral>
                FixedQuadrupleListTabulatedDihedral;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedDihedral::registerPython() {
            using namespace espresso::python;
            
            class_ <TabulatedDihedral, bases <DihedralPotential> >
                ("interaction_TabulatedDihedral", init <const char*>())
                .add_property("filename", &TabulatedDihedral::getFilename, &TabulatedDihedral::setFilename);
            
            class_ <FixedQuadrupleListTabulatedDihedral, bases <Interaction> > 
                ("interaction_FixedQuadrupleListTabulatedDihedral",
                init <shared_ptr<System>, shared_ptr <FixedQuadrupleList> >())
                .def("setPotential", &FixedQuadrupleListTabulatedDihedral::setPotential);
        }
        
    } // ns interaction
} // ns espresso