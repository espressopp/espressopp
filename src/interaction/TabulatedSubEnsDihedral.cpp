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
#include "TabulatedSubEnsDihedral.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"
#include "FixedQuadrupleListTypesInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {

        void TabulatedSubEnsDihedral::setFilenames(int dim,
            int itype, boost::python::list _filenames) {
            boost::mpi::communicator world;
            filenames.resize(dim);
            colVarRef.setDimension(dim);
            numInteractions = dim;
            for (int i=0; i<dim; ++i) {
              filenames[i] = boost::python::extract<std::string>(_filenames[i]);
              colVarRef[i].setDimension(8);
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

        void TabulatedSubEnsDihedral::addInteraction(int itype,
            boost::python::str fname, const RealND& _cvref) {
            boost::mpi::communicator world;
            int i = numInteractions;
            numInteractions += 1;
            colVarRef.setDimension(numInteractions);
            // Dimension 6: dihed,  bond, angle,
            // sd_dihed, sd_bond, sd_angle
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

        void TabulatedSubEnsDihedral::setColVarRef(
            const RealNDs& cvRefs){
            // Set the reference values of the collective variables
            // aka cluster centers
            for (int i=0; i<numInteractions; ++i)
                colVarRef[i] = cvRefs[i];
        }

        void TabulatedSubEnsDihedral::computeColVarWeights(
            const Real3D& dist21, const Real3D& dist32,
            const Real3D& dist43, const bc::BC& bc){
            // Sanity Check
            if (weights.getDimension() == 0)
                throw std::runtime_error("TabulatedSubEnsDihedral requires at least one interaction.");
            // Compute the weights for each force field
            // given the reference and instantaneous values of ColVars
            setColVar(dist21, dist32, dist43, bc);
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
                        // Choose between dihed, bonds, and angles
                        if (j <= 0+colVarDihedListSize) k = 0;
                        else if (j<1+colVarDihedListSize+colVarBondListSize) k = 1;
                        else k = 2;
                        real diff = colVar[j] -  colVarRef[i][k];
                        if (k == 0) {
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
        void TabulatedSubEnsDihedral::setColVar(
            const Real3D& r21, const Real3D& r32,
            const Real3D& r43, const bc::BC& bc) {
            colVar.setDimension(1);
            colVar[0] = computePhi(r21, r32, r43);
            int i=1;
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
            if (colVarBondList != nullptr) {
                // Now all bonds in colVarBondList
                colVarBondListSize = colVarBondList->size();
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
            if (colVarAngleList != nullptr) {
                // Now all angles in colVarAngleList
                colVarAngleListSize = colVarAngleList->size();
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
        }


        typedef class FixedQuadrupleListInteractionTemplate <TabulatedSubEnsDihedral>
            FixedQuadrupleListTabulatedSubEnsDihedral;

        typedef class FixedQuadrupleListTypesInteractionTemplate<TabulatedSubEnsDihedral>
            FixedQuadrupleListTypesTabulatedSubEnsDihedral;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedSubEnsDihedral::registerPython() {
            using namespace espressopp::python;

            class_ <TabulatedSubEnsDihedral, bases <DihedralPotential> >
                ("interaction_TabulatedSubEnsDihedral", init <>())
                .def("dimension_get", &TabulatedSubEnsDihedral::getDimension)
                .def("filenames_get", &TabulatedSubEnsDihedral::getFilenames)
                .def("filename_get", &TabulatedSubEnsDihedral::getFilename)
                .def("filename_set", &TabulatedSubEnsDihedral::setFilename)
                .def("targetProb_get", &TabulatedSubEnsDihedral::getTargetProb)
                .def("targetProb_set", &TabulatedSubEnsDihedral::setTargetProb)
                .def("colVarSd_get", &TabulatedSubEnsDihedral::getColVarSds)
                .def("colVarSd_set", &TabulatedSubEnsDihedral::setColVarSd)
                .def("weight_get", &TabulatedSubEnsDihedral::getWeights)
                .def("weight_set", &TabulatedSubEnsDihedral::setWeight)
                .def("alpha_get", &TabulatedSubEnsDihedral::getAlpha)
                .def("alpha_set", &TabulatedSubEnsDihedral::setAlpha)
                .def("addInteraction", &TabulatedSubEnsDihedral::addInteraction)
                .def("colVarRefs_get", &TabulatedSubEnsDihedral::getColVarRefs)
                .def("colVarRef_get", &TabulatedSubEnsDihedral::getColVarRef)
                .def_pickle(TabulatedSubEnsDihedral_pickle())
                ;

            class_ <FixedQuadrupleListTabulatedSubEnsDihedral, bases <Interaction> >
                ("interaction_FixedQuadrupleListTabulatedSubEnsDihedral",
                        init <shared_ptr<System>,
                              shared_ptr<FixedQuadrupleList>,
                              shared_ptr<TabulatedSubEnsDihedral> >())
                .def("setPotential", &FixedQuadrupleListTabulatedSubEnsDihedral::setPotential)
                .def("getFixedQuadrupleList", &FixedQuadrupleListTabulatedSubEnsDihedral::getFixedQuadrupleList);

            class_< FixedQuadrupleListTypesTabulatedSubEnsDihedral, bases< Interaction > >
                ("interaction_FixedQuadrupleListTypesTabulatedSubEnsDihedral",
                 init< shared_ptr<System>, shared_ptr<FixedQuadrupleList> >())
                .def("setPotential", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::setPotential)
                .def("getPotential", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::getPotentialPtr)
                .def("setFixedQuadrupleList", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::setFixedQuadrupleList)
                .def("getFixedQuadrupleList", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::getFixedQuadrupleList);

        }

    } // ns interaction
} // ns espressopp
