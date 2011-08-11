#include "python.hpp"
#include "TabulatedDihedral.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"

namespace espresso {
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

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedDihedral::registerPython() {
            using namespace espresso::python;
            
            class_ <TabulatedDihedral, bases <DihedralPotential> >
                ("interaction_TabulatedDihedral", init <int, const char*>())
                .add_property("filename", &TabulatedDihedral::getFilename, &TabulatedDihedral::setFilename);
            
            class_ <FixedQuadrupleListTabulatedDihedral, bases <Interaction> > 
                ("interaction_FixedQuadrupleListTabulatedDihedral",
                        init <shared_ptr<System>,
                              shared_ptr<FixedQuadrupleList>,
                              shared_ptr<TabulatedDihedral> >())
                .def("setPotential", &FixedQuadrupleListTabulatedDihedral::setPotential);
        }
        
    } // ns interaction
} // ns espresso
