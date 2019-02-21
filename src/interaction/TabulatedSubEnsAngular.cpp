/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research

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
#include "DihedralPotential.hpp"

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
            colVarRef.setDimension(numInteractions);
            // Dimension 6: angle, bond, dihed,
            // sd_angle, sd_bond, sd_dihed
            colVarRef[i].setDimension(6);
            colVarRef[i] = _cvref;
            filenames.push_back(boost::python::extract<std::string>(fname));
            weights.push_back(0.);
            weightSum.push_back(0.);
            targetProb.push_back(0.);
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

        void TabulatedSubEnsAngular::setColVarRef(
            const RealNDs& cvRefs){
            // Set the reference values of the collective variables
            // aka cluster centers
            for (int i=0; i<numInteractions; ++i)
                colVarRef[i] = cvRefs[i];
        }

        void TabulatedSubEnsAngular::computeColVarWeights(
            const Real3D& dist12, const Real3D& dist32, const bc::BC& bc){
            // Sanity Check
            if (weights.getDimension() == 0)
                throw std::runtime_error("TabulatedSubEnsAngular requires at least one interaction.");
            // Compute the weights for each force field
            // given the reference and instantaneous values of ColVars
            setColVar(dist12, dist32, bc);
            // Compute weights up to next to last FF
            real maxWeight = 0.;
            int maxWeightI = 0;
            // Check first whether we're stuck in a surface
            bool stuck = false;
            for (int i=0; i<numInteractions; ++i) {
                if (weights[i] > maxWeight) {
                    maxWeight = weights[i];
                    maxWeightI = i;
                }
            }
            // Lock the surface until it's sampled >= 98% of the target prob.
            if (weightCounts > 0 &&
                maxWeightI < numInteractions-1 &&
                maxWeight > 0. &&
                weightSum[maxWeightI]/weightCounts < 0.98*targetProb[maxWeightI])
                stuck = true;
            if (!stuck) {
                maxWeight = 0.;
                maxWeightI = numInteractions-1;
                for (int i=0; i<numInteractions-1; ++i) {
                    weights[i]    = 1.;
                    real norm_d_i = 0.;
                    real norm_l_i = 0.;
                    for (int j=0; j<colVar.getDimension(); ++j) {
                        int k = 0;
                        // Choose between angle, bond, and dihed
                        if (j <= 0+colVarAngleListSize) k = 0;
                        else if (j<1+colVarAngleListSize+colVarBondListSize) k = 1;
                        else k = 2;
                        real diff = colVar[j] -  colVarRef[i][k];
                        if (k == 2) {
                            if (diff>M_PI) diff -= 2.0*M_PI;
                            if (diff<(-1.0*M_PI)) diff += 2.0*M_PI;
                        }
                        norm_d_i += pow(diff / colVarSd[k], 2);
                        norm_l_i += pow(colVarRef[i][3+k], 2);
                    }
                    if (norm_d_i > norm_l_i)
                      weights[i] = exp(- (sqrt(norm_d_i) - sqrt(norm_l_i)) / alpha);
                    if (weights[i] > maxWeight) {
                      maxWeight = weights[i];
                      maxWeightI = i;
                    }
                }
                for (int i=0; i<numInteractions-1; ++i) {
                    if (i != maxWeightI)
                        weights[i] = 0.;
                    else {
                        if (weightCounts > 0 &&
                            weights[i] > 0.01 &&
                            weightSum[i]/weightCounts < 0.98*targetProb[i]) {
                            weights[i] = 1.;
                            maxWeight = 1.;
                        }
                    }
                }
                if (maxWeightI == numInteractions-1)
                    maxWeight = 1.;
                weights[numInteractions-1] = 1. - maxWeight;
            }

            // Update weightSum
            for (int i=0; i<numInteractions; ++i)
                weightSum[i] += weights[i];
            weightCounts += 1;
        }

        // Collective variables
        void TabulatedSubEnsAngular::setColVar(const Real3D& dist12,
              const Real3D& dist32, const bc::BC& bc) {
            colVar.setDimension(1);
            real dist12_sqr = dist12 * dist12;
            real dist32_sqr = dist32 * dist32;
            real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
            real cos_theta = dist12 * dist32 / dist1232;
            colVar[0] = acos(cos_theta);
            int i=1;
            if (colVarAngleList != nullptr) {
                colVarAngleListSize = colVarAngleList->size();
                // Now all angles in colVarAngleList
                for (FixedTripleList::TripleList::Iterator it(*colVarAngleList); it.isValid(); ++it) {
                  colVar.setDimension(i+1);
                  Particle &p1 = *it->first;
                  Particle &p2 = *it->second;
                  Particle &p3 = *it->third;
                  Real3D dist12, dist32;
                  bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
                  bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
                  real dist12_sqr = dist12 * dist12;
                  real dist32_sqr = dist32 * dist32;
                  real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
                  real cos_theta = dist12 * dist32 / dist1232;
                  colVar[i] = acos(cos_theta);
                  i+=1;
                }
            }
            if (colVarBondList != nullptr) {
                colVarBondListSize = colVarBondList->size();
                // Now all bonds in colVarBondList
                for (FixedPairList::PairList::Iterator it(*colVarBondList); it.isValid(); ++it) {
                  colVar.setDimension(i+1);
                  Particle &p1 = *it->first;
                  Particle &p2 = *it->second;
                  Real3D dist12;
                  bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
                  colVar[i] = sqrt(dist12 * dist12);
                  i+=1;
                }
            }
            if (colVarDihedList != nullptr) {
                colVarDihedListSize = colVarDihedList->size();
                // Now all dihedrals in colVarDihedList
                for (FixedQuadrupleList::QuadrupleList::Iterator it(*colVarDihedList); it.isValid(); ++it) {
                  colVar.setDimension(i+1);
                  Particle &p1 = *it->first;
                  Particle &p2 = *it->second;
                  Particle &p3 = *it->third;
                  Particle &p4 = *it->fourth;
                  Real3D r21, r32, r43;
                  bc.getMinimumImageVectorBox(r21, p2.position(), p1.position());
                  bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
                  bc.getMinimumImageVectorBox(r43, p4.position(), p3.position());
                  colVar[i] = DihedralPotential::computePhi(r21, r32, r43);
                  i+=1;
                }
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
                ("interaction_TabulatedSubEnsAngular", init <>())
                .def("dimension_get", &TabulatedSubEnsAngular::getDimension)
                .def("filenames_get", &TabulatedSubEnsAngular::getFilenames)
                .def("filename_get", &TabulatedSubEnsAngular::getFilename)
                .def("filename_set", &TabulatedSubEnsAngular::setFilename)
                .def("targetProb_get", &TabulatedSubEnsAngular::getTargetProb)
                .def("targetProb_set", &TabulatedSubEnsAngular::setTargetProb)
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
