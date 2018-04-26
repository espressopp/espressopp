/*
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
#include "TabulatedSubEnsAngular.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedTripleListTypesInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {

        void TabulatedSubEnsAngular::setFilename(int itype, char** _filenames) {
            boost::mpi::communicator world;
            for (int i=0; i<numInteractions; ++i) {
              if (itype == 1) { // create a new InterpolationLinear
                  tables[i] = make_shared <InterpolationLinear> ();
                  tables[i]->read(world, _filenames[i]);
              }

              else if (itype == 2) { // create a new InterpolationAkima
                  tables[i] = make_shared <InterpolationAkima> ();
                  tables[i]->read(world, _filenames[i]);
              }

              else if (itype == 3) { // create a new InterpolationCubic
                  tables[i] = make_shared <InterpolationCubic> ();
                  tables[i]->read(world, _filenames[i]);
              }
          }
        }

        double TabulatedSubEnsAngular::distColVars(
            std::array<real, 4> cv1, std::array<real, 4> cv2){
            // Compute distance between colvars cv1 and cv2
            // Metric is euclidean distance
            real dist = 0.;
            for (int i=0; i<3; ++i)
                dist += pow(cv1[i] - cv2[i], 2);
            return sqrt(dist);
        }

        void TabulatedSubEnsAngular::setColVarRef(
            std::vector<std::array<real, 4>> cvRefs){
            // Set the reference values of the collective variables
            // aka cluster centers
            for (int i=0; i<numInteractions; ++i) {
                for (int j=0; j<4; ++j)
                    colVarRef[i][j] = cvRefs[i][j];
                }
        }

        void TabulatedSubEnsAngular::computeWeights(){
            // Compute the weights for each force field
            // given the reference and instantaneous values of ColVars
            real normalization = 0.;
            for (int i=0; i<numInteractions; ++i) {
                // Distance between inst and ref_i
                setColVarRef(colVarRef);
                real d_i = distColVars(colVarInst, colVarRef[i]);
                // Lengthscale of the cluster i
                real length_ci = colVarRef[i][3];
                if (d_i < length_ci) weights[i] = 1.;
                else weights[i] = exp(- d_i / (alpha * length_ci));
                normalization += weights[i];
            }
            // Normalize
            for (int i=0; i<numInteractions; ++i) {
                if (weights[i] > 0.)
                    weights[i] /= normalization;
            }
        }

        typedef class FixedTripleListInteractionTemplate <TabulatedSubEnsAngular>
                FixedTripleListTabulatedSubEnsAngular;
        typedef class FixedTripleListTypesInteractionTemplate<TabulatedSubEnsAngular>
            FixedTripleListTypesTabulatedSubEnsAngular;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedSubEnsAngular::registerPython() {
            using namespace espressopp::python; 

            class_ <TabulatedSubEnsAngular, bases <AngularPotential> >
                ("interaction_TabulatedSubEnsAngular", init <int, int, char**>())
                .add_property("filename", &TabulatedSubEnsAngular::getFilename, &TabulatedSubEnsAngular::setFilename)
                .def_pickle(TabulatedSubEnsAngular_pickle());

            class_ <FixedTripleListTabulatedSubEnsAngular, bases <Interaction> >
                ("interaction_FixedTripleListTabulatedSubEnsAngular",
                init <shared_ptr<System>,
                      shared_ptr <FixedTripleList>,
                      shared_ptr <TabulatedSubEnsAngular> >())
                .def("setPotential", &FixedTripleListTabulatedSubEnsAngular::setPotential)
                .def("getFixedTripleList", &FixedTripleListTabulatedSubEnsAngular::getFixedTripleList);

            class_< FixedTripleListTypesTabulatedSubEnsAngular, bases< Interaction > >
                ("interaction_FixedTripleListTypesTabulatedSubEnsAngular",
                 init< shared_ptr<System>, shared_ptr<FixedTripleList> >())
                .def("setPotential", &FixedTripleListTypesTabulatedSubEnsAngular::setPotential)
                .def("getPotential", &FixedTripleListTypesTabulatedSubEnsAngular::getPotentialPtr)
                .def("setFixedTripleList", &FixedTripleListTypesTabulatedSubEnsAngular::setFixedTripleList)
                .def("getFixedTripleList", &FixedTripleListTypesTabulatedSubEnsAngular::getFixedTripleList);
        }

    } // ns interaction
} // ns espressopp
