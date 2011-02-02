#include "python.hpp"
#include "TabulatedAngular.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
    namespace interaction {
        
        void TabulatedAngular::setFilename(const char* _filename) {
            boost::mpi::communicator world;
            filename = _filename;
            // create a new InterpolationTable
            table = make_shared <InterpolationTable> ();
            table->read(world, _filename);
        }

        typedef class FixedTripleListInteractionTemplate <TabulatedAngular>
                FixedTripleListTabulatedAngular;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedAngular::registerPython() {
            using namespace espresso::python;
            
            class_ <TabulatedAngular, bases <AngularPotential> >
                ("interaction_TabulatedAngular", init <const char*>())
                .add_property("filename", &TabulatedAngular::getFilename, &TabulatedAngular::setFilename);
            
            class_ <FixedTripleListTabulatedAngular, bases <Interaction> > 
                ("interaction_FixedTripleListTabulatedAngular",
                init <shared_ptr<System>, shared_ptr <FixedTripleList> >())
                .def("setPotential", &FixedTripleListTabulatedAngular::setPotential);
        }
        
    } // ns interaction
} // ns espresso