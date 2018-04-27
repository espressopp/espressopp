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

        void TabulatedSubEnsAngular::setFilenames(int dim,
            int itype, std::vector<std::string> _filenames) {
            boost::mpi::communicator world;
            filenames = _filenames;
            numInteractions = dim;
            for (int i=0; i<dim; ++i) {
              if (itype == 1) { // create a new InterpolationLinear
                  tables[i] = make_shared <InterpolationLinear> ();
                  tables[i]->read(world, _filenames[i].c_str());
              }

              else if (itype == 2) { // create a new InterpolationAkima
                  tables[i] = make_shared <InterpolationAkima> ();
                  tables[i]->read(world, _filenames[i].c_str());
              }

              else if (itype == 3) { // create a new InterpolationCubic
                  tables[i] = make_shared <InterpolationCubic> ();
                  tables[i]->read(world, _filenames[i].c_str());
              }
          }
        }

        double TabulatedSubEnsAngular::distColVars(
            cv_angle cv1, cv_angle cv2){
            // Compute distance between colvars cv1 and cv2
            // Metric is euclidean distance
            real dist = 0.;
            dist = pow(cv1.b1 - cv2.b1, 2)
                 + pow(cv1.b2 - cv2.b2, 2)
                 + pow(cv1.theta - cv2.theta, 2);
            return sqrt(dist);
        }

        void TabulatedSubEnsAngular::setColVarRef(
            std::vector< cv_angle > cvRefs){
            // Set the reference values of the collective variables
            // aka cluster centers
            for (int i=0; i<numInteractions; ++i)
                colVarRef[i] = cvRefs[i];
        }

        void TabulatedSubEnsAngular::computeColVarWeights(
            const Real3D& dist12, const Real3D& dist32){
            // Compute the weights for each force field
            // given the reference and instantaneous values of ColVars
            setColVar(dist12, dist32);
            real normalization = 0.;
            for (int i=0; i<numInteractions; ++i) {
                // Distance between inst and ref_i
                real d_i = distColVars(colVar, colVarRef[i]);
                // Lengthscale of the cluster i
                real length_ci = colVarRef[i].size;
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
                ("interaction_TabulatedSubEnsAngular", init <int, int, std::vector<std::string>>())
                .add_property("numInteractions", &TabulatedSubEnsAngular::getDimension, &TabulatedSubEnsAngular::setDimension)
                .add_property("filenames", &TabulatedSubEnsAngular::getFilenames, &TabulatedSubEnsAngular::setFilenames)
                .add_property("colVarRef", &TabulatedSubEnsAngular::getColVarRef, &TabulatedSubEnsAngular::setColVarRef)
                .add_property("weights", &TabulatedSubEnsAngular::getWeights)
                .add_property("alpha", &TabulatedSubEnsAngular::getWeightScalingFactor, &TabulatedSubEnsAngular::setWeightScalingFactor)
                .def_pickle(TabulatedSubEnsAngular_pickle())
                ;

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
