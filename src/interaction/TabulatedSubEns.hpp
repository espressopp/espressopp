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

// ESPP_CLASS
#ifndef _INTERACTION_TABULATEDSUBENS_HPP
#define _INTERACTION_TABULATEDSUBENS_HPP

//#include <stdexcept>
#include "Potential.hpp"
#include "Interpolation.hpp"
#include "RealND.hpp"
#include "bc/BC.hpp"

namespace espressopp {

  namespace interaction {

    /** This class provides methods to compute forces and energies of
	a tabulated potential.

        The potential and forces must be provided in a file.

        Be careful: default and copy constructor of this class are used.
    */

    class TabulatedSubEns: public PotentialTemplate <TabulatedSubEns> {

        private:
            int numInteractions;
            std::vector<std::string> filenames;
            std::vector<shared_ptr <Interpolation>> tables;
            int interpolationType;
            // Reference values of the collective variable centers
            RealNDs colVarRef;
            // Weights of each table
            RealND weights;
            // Target probability of each table
            RealND targetProb;
            // Running sum of each weight and number of counts
            RealND weightSum;
            int weightCounts;
            // Renormalize collective variables: std
            RealND colVarSd;
            // characteristic decay length of the interpolation
            real alpha;
            // Size of CV partners
            int colVarBondListSize;
            int colVarAngleListSize;
            int colVarDihedListSize;

        public:
            static void registerPython();

            TabulatedSubEns() :
              numInteractions(0) {
              setCutoff(infinity);
              weights.setDimension(0);
              weightSum.setDimension(0);
              targetProb.setDimension(0);
              weightCounts = 0;
              colVarSd.setDimension(3);
              colVarRef.setDimension(0);
              alpha = 1.;
              colVarBondListSize = 0;
              colVarAngleListSize = 0;
              colVarDihedListSize = 0;
            }

            void addInteraction(int itype, boost::python::str fname,
                                const RealND& _cvref);

            void setDimension(int _dim) {
              numInteractions = _dim;
              colVarRef.setDimension( numInteractions );
              tables.resize( numInteractions );
              filenames.resize( numInteractions );
              weights.setDimension( numInteractions );
              weightSum.setDimension( numInteractions );
              targetProb.setDimension( numInteractions );
            }

            int getDimension() const { return numInteractions; }

            /** Setter for the interpolation type */
            void setInterpolationType(int itype) { interpolationType = itype; }

            /** Getter for the interpolation type */
            int getInterpolationType() const { return interpolationType; }

            RealND getTargetProb() const { return targetProb; }

            void setTargetProb(int index, real _r) {
                return targetProb.setItem(index, _r);
            }

            RealND getColVarSds() const { return colVarSd; }

            void setColVarSd(int index, real _r) {
                return colVarSd.setItem(index, _r);
            }

            RealND getColVarRef(int i) const { return colVarRef[i]; }

            RealNDs getColVarRefs() const { return colVarRef; }

            void setColVarRef(const RealNDs& cvRefs);

            void setColVarRefs(const RealNDs& c) { colVarRef = c; }

            boost::python::list getFilenames() const {
                return boost::python::list(filenames); }

            boost::python::list getFilename(int index) const {
                return boost::python::list(filenames[index]);
            }

            void setFilename(int index, boost::python::list _f) {
                filenames[index] = boost::python::extract<std::string>(_f);
            }

            void setFilenames(int dim, int itype, boost::python::list _filenames);

            RealND getWeights() const { return weights; }

            void setWeights(const RealND& r) { weights = r; }

            real getWeight(int index) const {
                return weights.getItem(index);
            }

            void setWeight(int index, real _w) {
                return weights.setItem(index, _w);
            }

            real getAlpha() const { return alpha; }

            void setAlpha(real _r) { alpha = _r; }

            void computeColVarWeights(const Real3D& dist, const bc::BC& bc);

            void setColVar(const Real3D& dist, const bc::BC& bc);

            real _computeEnergySqrRaw(real distSqr) const {
              real e = 0.;
              for (int i=0; i<numInteractions; ++i)
                  e += weights[i] * tables[i]->getEnergy(sqrt(distSqr));
              return e;
            }

            bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
                real ffactor = 0;
                real distrt = sqrt(distSqr);
                for	(int i=0; i<numInteractions; ++i)
                    ffactor += weights[i] * tables[i]->getForce(distrt) / distrt;
                force = dist * ffactor;
                return true;
            }

    };//class

    // provide pickle support
    struct TabulatedSubEns_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(TabulatedSubEns const& pot)
      {
        int itp = pot.getInterpolationType();
        boost::python::list fns;
        RealNDs cvrefs = pot.getColVarRefs();
        int dim = pot.getDimension();
        fns = pot.getFilenames();
        RealND cvsd = pot.getColVarSds();
        real rc = pot.getCutoff();
        real alp = pot.getAlpha();
        return boost::python::make_tuple(dim, itp, fns, cvrefs,
                                         cvsd, alp, rc);
      }
    };

  }
}

#endif
