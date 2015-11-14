/*
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
#include "TabulatedAngular.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedTripleListPIadressInteractionTemplate.hpp"

namespace espressopp {
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
        typedef class FixedTripleListPIadressInteractionTemplate <TabulatedAngular>
                FixedTripleListPIadressTabulatedAngular;
        
        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedAngular::registerPython() {
            using namespace espressopp::python;
            
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
   
            class_ <FixedTripleListPIadressTabulatedAngular, bases <Interaction> > 
                ("interaction_FixedTripleListPIadressTabulatedAngular",
                init <shared_ptr<System>,
                      shared_ptr <FixedTripleList>,
                      shared_ptr<FixedTupleListAdress>,
                      shared_ptr <TabulatedAngular>,
                      int>())
                .def("setPotential", &FixedTripleListPIadressTabulatedAngular::setPotential)
                .def("getFixedTripleList", &FixedTripleListPIadressTabulatedAngular::getFixedTripleList);            
        }
        
    } // ns interaction
} // ns espressopp
