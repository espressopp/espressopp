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
#include "TabulatedSubEns.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {


    void TabulatedSubEns::setFilenames(int dim,
      int itype, boost::python::list _filenames) {
        boost::mpi::communicator world;
        filenames.resize(dim);
        colVarRef.setDimension(dim);
        numInteractions = dim;
        for (int i=0; i<dim; ++i) {
          filenames[i] = boost::python::extract<std::string>(_filenames[i]);
          colVarRef[i].setDimension(4);
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

    void TabulatedSubEns::addInteraction(int itype,
        boost::python::str fname, const RealND& _cvref) {
        boost::mpi::communicator world;
        int i = numInteractions;
        numInteractions += 1;
        colVarRef.setDimension(numInteractions);
        // Dimension 6: angle, bond, dihed, sd_angle, sd_bond, sd_dihed
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

    void TabulatedSubEns::setColVarRef(
        const RealNDs& cvRefs){
        // Set the reference values of the collective variables
        // aka cluster centers
        for (int i=0; i<numInteractions; ++i)
            colVarRef[i] = cvRefs[i];
    }

    void TabulatedSubEns::computeColVarWeights(const Real3D& dist,
        const bc::BC& bc){
        // Compute the weights for each force field
        // given the reference and instantaneous values of ColVars
        setColVar(dist, bc);
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
        if (weightCounts > 0 &&
            maxWeightI < numInteractions-1 &&
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
                    // Choose between bond, angle, and dihed
                    if (j <= 0+colVarBondList->size()) k = 0;
                    else if (j>0 && j<1+colVarBondList->size()+colVarAngleList->size()) k = 1;
                    else k = 2;
                    norm_d_i += pow((colVar[j] -  colVarRef[i][k]) / colVarSd[k], 2);
                    norm_l_i += pow(colVarRef[i][3+k], 2);
                }
                if (norm_d_i > norm_l_i)
                  weights[i]  = exp(- (sqrt(norm_d_i) - sqrt(norm_l_i)) / alpha);
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
                        weights[i] = 1.0;
                        maxWeight = weights[i];
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
    void TabulatedSubEns::setColVar(const Real3D& dist, const bc::BC& bc) {
        colVar.setDimension(1+colVarBondList->size()+colVarAngleList->size());
        colVar[0] = sqrt(dist*dist);
        // Now all bonds in colVarBondList
        int i=1;
        for (FixedPairList::PairList::Iterator it(*colVarBondList); it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          Real3D dist12;
          bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
          colVar[i] = sqrt(dist12 * dist12);
          i+=1;
        }
        // Now all angles in colVarAngleList
        for (FixedTripleList::TripleList::Iterator it(*colVarAngleList); it.isValid(); ++it) {
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

    typedef class VerletListInteractionTemplate <TabulatedSubEns> VerletListTabulatedSubEns;
    typedef class VerletListAdressInteractionTemplate <TabulatedSubEns, TabulatedSubEns> VerletListAdressTabulatedSubEns;
    typedef class VerletListHadressInteractionTemplate <TabulatedSubEns, TabulatedSubEns> VerletListHadressTabulatedSubEns;
    typedef class CellListAllPairsInteractionTemplate <TabulatedSubEns> CellListTabulatedSubEns;
    typedef class FixedPairListInteractionTemplate <TabulatedSubEns> FixedPairListTabulatedSubEns;
    typedef class FixedPairListTypesInteractionTemplate <TabulatedSubEns> FixedPairListTypesTabulatedSubEns;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void TabulatedSubEns::registerPython() {
      using namespace espressopp::python;

      class_ <TabulatedSubEns, bases <Potential> >
        ("interaction_TabulatedSubEns", init <>())
            .def("dimension_get", &TabulatedSubEns::getDimension)
            .def("filenames_get", &TabulatedSubEns::getFilenames)
            .def("filename_get", &TabulatedSubEns::getFilename)
            .def("filename_set", &TabulatedSubEns::setFilename)
            .def("targetProb_get", &TabulatedSubEns::getTargetProb)
            .def("targetProb_set", &TabulatedSubEns::setTargetProb)
            .def("colVarMu_get", &TabulatedSubEns::getColVarMus)
            .def("colVarMu_set", &TabulatedSubEns::setColVarMu)
            .def("colVarSd_get", &TabulatedSubEns::getColVarSds)
            .def("colVarSd_set", &TabulatedSubEns::setColVarSd)
            .def("weight_get", &TabulatedSubEns::getWeights)
            .def("weight_set", &TabulatedSubEns::setWeight)
            .def("alpha_get", &TabulatedSubEns::getAlpha)
            .def("alpha_set", &TabulatedSubEns::setAlpha)
            .def("addInteraction", &TabulatedSubEns::addInteraction)
            .def("colVarRefs_get", &TabulatedSubEns::getColVarRefs)
            .def("colVarRef_get", &TabulatedSubEns::getColVarRef)
            .def_pickle(TabulatedSubEns_pickle())
        ;

      class_ <VerletListTabulatedSubEns, bases <Interaction> >
        ("interaction_VerletListTabulatedSubEns", init <shared_ptr<VerletList> >())
            .def("setPotential", &VerletListTabulatedSubEns::setPotential)
            .def("getPotential", &VerletListTabulatedSubEns::getPotentialPtr)
        ;

      class_ <VerletListAdressTabulatedSubEns, bases <Interaction> >
        ("interaction_VerletListAdressTabulatedSubEns",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress> >()
                )
            .def("setPotentialAT", &VerletListAdressTabulatedSubEns::setPotentialAT)
            .def("setPotentialCG", &VerletListAdressTabulatedSubEns::setPotentialCG);
        ;

      class_ <VerletListHadressTabulatedSubEns, bases <Interaction> >
        ("interaction_VerletListHadressTabulatedSubEns",
           init <shared_ptr<VerletListAdress>,
                 shared_ptr<FixedTupleListAdress> >()
                )
            .def("setPotentialAT", &VerletListHadressTabulatedSubEns::setPotentialAT)
            .def("setPotentialCG", &VerletListHadressTabulatedSubEns::setPotentialCG);
        ;

      class_ <CellListTabulatedSubEns, bases <Interaction> >
        ("interaction_CellListTabulatedSubEns", init <shared_ptr <storage::Storage> >())
            .def("setPotential", &CellListTabulatedSubEns::setPotential);
        ;

      class_ <FixedPairListTabulatedSubEns, bases <Interaction> >
        ("interaction_FixedPairListTabulatedSubEns",
          init <shared_ptr<System>,
                shared_ptr<FixedPairList>,
                shared_ptr<TabulatedSubEns> >()
        )
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<TabulatedSubEns> >())
        .def("setPotential", &FixedPairListTabulatedSubEns::setPotential)
        .def("setFixedPairList", &FixedPairListTabulatedSubEns::setFixedPairList)
        .def("getFixedPairList", &FixedPairListTabulatedSubEns::getFixedPairList);

      class_< FixedPairListTypesTabulatedSubEns, bases< Interaction > >
          ("interaction_FixedPairListTypesTabulatedSubEns",
           init< shared_ptr<System>, shared_ptr<FixedPairList> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
          .def("setPotential", &FixedPairListTypesTabulatedSubEns::setPotential)
          .def("getPotential", &FixedPairListTypesTabulatedSubEns::getPotentialPtr)
          .def("setFixedPairList", &FixedPairListTypesTabulatedSubEns::setFixedPairList)
          .def("getFixedPairList", &FixedPairListTypesTabulatedSubEns::getFixedPairList);
    }

  }
}
