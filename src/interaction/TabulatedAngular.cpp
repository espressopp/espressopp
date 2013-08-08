#include "python.hpp"
#include "TabulatedAngular.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
    namespace interaction {
        
        void TabulatedAngular::setFilename(int itype, const char* _filename) {
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

        typedef class FixedTripleListInteractionTemplate <TabulatedAngular>
                FixedTripleListTabulatedAngular;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedAngular::registerPython() {
            using namespace espresso::python;
            
            class_ <TabulatedAngular, bases <AngularPotential> >
                ("interaction_TabulatedAngular", init <int, const char*>())
                .add_property("filename", &TabulatedAngular::getFilename, &TabulatedAngular::setFilename);
            
            class_ <FixedTripleListTabulatedAngular, bases <Interaction> > 
                ("interaction_FixedTripleListTabulatedAngular",
                init <shared_ptr<System>,
                      shared_ptr <FixedTripleList>,
                      shared_ptr <TabulatedAngular> >())
                .def("setPotential", &FixedTripleListTabulatedAngular::setPotential)
                .def("getFixedTripleList", &FixedTripleListTabulatedAngular::getFixedTripleList);
        }
        
    } // ns interaction
} // ns espresso
