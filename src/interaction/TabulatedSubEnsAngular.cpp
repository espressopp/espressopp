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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace espressopp {
    namespace interaction {

        void TabulatedSubEnsAngular::setFilenames(int dim,
            int itype, boost::python::list _filenames) {
            boost::mpi::communicator world;
            filenames.resize(dim);
            colVarRef.setDimension(dim);
            numInteractions = dim;
            for (int i=0; i<dim; ++i) {
              filenames[i] = boost::python::extract<std::string>(_filenames[i]);
              colVarRef[i].setDimension(6);
              if (itype == 1) { // create a new InterpolationLinear
                  tables[i] = make_shared <InterpolationLinear> ();
                  tables[i]->read(world, filenames[i].c_str());
              }

              else if (itype == 2) { // create a new InterpolationAkima
                  tables[i] = make_shared <InterpolationAkima> ();
                  tables[i]->read(world, filenames[i].c_str());
              }

              else if (itype == 3) { // create a new InterpolationCubic
                  tables[i] = make_shared <InterpolationCubic> ();
                  tables[i]->read(world, filenames[i].c_str());
              }
          }
        }

        void TabulatedSubEnsAngular::addInteraction(int itype,
            boost::python::str fname, const RealND& _cvref) {
            boost::mpi::communicator world;
            int i = numInteractions;
            numInteractions += 1;
            colVarRef.push_back(_cvref);
            filenames.push_back(boost::python::extract<std::string>(fname));
            weights.push_back(0.);
            if (itype == 1) { // create a new InterpolationLinear
                  tables.push_back(make_shared <InterpolationLinear> ());
                  tables[i]->read(world, filenames[i].c_str());
              }
              else if (itype == 2) { // create a new InterpolationAkima
                  tables.push_back(make_shared <InterpolationAkima> ());
                  tables[i]->read(world, filenames[i].c_str());
              }
              else if (itype == 3) { // create a new InterpolationCubic
                  tables.push_back(make_shared <InterpolationCubic> ());
                  tables[i]->read(world, filenames[i].c_str());
              }
        }

        double TabulatedSubEnsAngular::distColVars(
            const RealND& cv1, const RealND& cv2, int i){
            // Compute distance between colvars cv1 and cv2
            // in renormalized space
            // Metric is euclidean distance
            real dist = 0.;
            // Try with only one bond and one angle
            // for (int i = 0; i<3; ++i)
            //     dist += pow((cv1[i] - cv2[i]) / colVarSd[i], 2);
            // return sqrt(dist);
            return abs((cv1[i] - cv2[i]) / colVarSd[i]);
        }

        void TabulatedSubEnsAngular::setColVarRef(
            const RealNDs& cvRefs){
            // Set the reference values of the collective variables
            // aka cluster centers
            for (int i=0; i<numInteractions; ++i)
                colVarRef[i] = cvRefs[i];
        }

        void TabulatedSubEnsAngular::computeColVarWeights(
            const Real3D& dist12, const Real3D& dist32){
            // Compute the weights for each force field
            // given the reference and instantaneous values of ColVars
            std::ostringstream msg;
            setColVar(dist12, dist32);
            real sumWeights = 0.;
            // Compute weights up to next to last FF
            for (int i=0; i<numInteractions-1; ++i) {
                weights[i] = 1.;
                for (int j=0; j<3; ++j) {
                    // Lengthscale of the cluster i
                    real length_ci = colVarRef[i][3+j];
                    // Distance between inst and ref_i
                    real d_i = distColVars(colVar, colVarRef[i], j);
                    if (d_i > length_ci)
                        weights[i] *= exp(- (d_i - length_ci) / alpha);
                }

                real d_i = distColVars(colVar, colVarRef[i], 0);
                real length_ci = colVarRef[i][3];
                // if (d_i < length_ci) 
                //   weights[i] = 1.;
                // else weights[i] = exp(- (d_i - length_ci) / ( alpha ));
                // std::cout << "length_ci " << length_ci << " weight " << i
                //     << " " << weights[i] << " di " << d_i << std::endl;
                sumWeights += weights[i];
            }
            // real sumWeightsNew = sumWeights;
            // if (sumWeights > 1.) {
            //   // Rescale down
            //   sumWeightsNew = 0.;
            //   for (int j=0; j<numInteractions-1; ++j) {
            //     weights[j] = std::max(0., weights[j]-(sumWeights-1.)/(numInteractions-1));
            //     sumWeightsNew += weights[j];
            //   }
            // }
            // Last FF gets the rest of the contribution
            // weights[numInteractions-1] = 1. - sumWeightsNew;
            weights[numInteractions-1] = 1.;
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

            class_<VectorStrings>("VectorStrings")
                .def(vector_indexing_suite<VectorStrings>() );

            class_ <TabulatedSubEnsAngular, bases <AngularPotential> >
                ("interaction_TabulatedSubEnsAngular", init <>())
                .def("dimension_get", &TabulatedSubEnsAngular::getDimension)
                .def("filenames_get", &TabulatedSubEnsAngular::getFilenames)
                .def("filename_get", &TabulatedSubEnsAngular::getFilename)
                .def("filename_set", &TabulatedSubEnsAngular::setFilename)
                .def("colVarMu_get", &TabulatedSubEnsAngular::getColVarMus)
                .def("colVarMu_set", &TabulatedSubEnsAngular::setColVarMu)
                .def("colVarSd_get", &TabulatedSubEnsAngular::getColVarSds)
                .def("colVarSd_set", &TabulatedSubEnsAngular::setColVarSd)
                .def("weight_get", &TabulatedSubEnsAngular::getWeights)
                .def("weight_set", &TabulatedSubEnsAngular::setWeight)
                .def("alpha_get", &TabulatedSubEnsAngular::getAlpha)
                .def("alpha_set", &TabulatedSubEnsAngular::setAlpha)
                .def("addInteraction", &TabulatedSubEnsAngular::addInteraction)
                .def("colVarRefs_get", &TabulatedSubEnsAngular::getColVarRefs)
                .def("colVarRef_get", &TabulatedSubEnsAngular::getColVarRef)
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
