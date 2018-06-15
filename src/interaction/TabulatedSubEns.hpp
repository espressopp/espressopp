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
            // Weights and derivatives of each table
            RealND weights;
            RealND dweights;
            // Renormalize collective variables: mean and std
            RealND colVarMu;
            RealND colVarSd;
            // characteristic decay length of the interpolation
            real alpha;
            // offset free energy of each table
            RealND offsets;

        public:
            static void registerPython();

            TabulatedSubEns() :
              numInteractions(0) {
              setCutoff(infinity);
              weights.setDimension(0);
              dweights.setDimension(0);
              colVarMu.setDimension(3);
              colVarSd.setDimension(3);
              colVarRef.setDimension(0);
              alpha = 1.;
              offsets.setDimension(0);
            }

            void addInteraction(int itype, boost::python::str fname,
                                const RealND& _cvref, real _offset);

            void setDimension(int _dim) {
              numInteractions = _dim;
              colVarRef.setDimension( numInteractions );
              tables.resize( numInteractions );
              filenames.resize( numInteractions );
              weights.setDimension( numInteractions );
              dweights.setDimension( numInteractions );
              offsets.setDimension( numInteractions );
            }

            int getDimension() const { return numInteractions; }

            /** Setter for the interpolation type */
            void setInterpolationType(int itype) { interpolationType = itype; }

            /** Getter for the interpolation type */
            int getInterpolationType() const { return interpolationType; }

            RealND getColVarMus() const { return colVarMu; }

            void setColVarMu(int index, real _r) {
                return colVarMu.setItem(index, _r);
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

            RealND getDWeights() const { return dweights; }

            void setDWeight(int index, real _dw) {
                return dweights.setItem(index, _dw);
            }

            RealND getOffsets() const { return offsets; }

            void setOffset(int index, real _o) {
                return offsets.setItem(index, _o);
            }

            real getAlpha() const { return alpha; }

            void setAlpha(real _r) { alpha = _r; }

            void computeColVarWeights(const Real3D& dist, const bc::BC& bc);

            void setColVar(const Real3D& dist, const bc::BC& bc);

            real _computeEnergySqrRaw(real distSqr) const {
              real e = 0.;
              for	(int i=0; i<numInteractions; ++i) {
                  e += weights[i] * (tables[i]->getEnergy(sqrt(distSqr))
                          + offsets[i]);
              }
              return e;
            }

            bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
                real ffactor = 0;
                real distrt = sqrt(distSqr);
                for	(int i=0; i<numInteractions; ++i)
                    ffactor += weights[i] * tables[i]->getForce(distrt) / distrt
                              - offsets[i] * dweights[i] / distrt;
                force = dist * ffactor;
                return true;
            }

            real computeForceNorm(real d, RealND w, RealND dw) {
              weights = w;
              dweights = dw;
              real ffactor = 0;
              real distrt = d;
              for	(int i=0; i<numInteractions; ++i)
                  ffactor += weights[i] * tables[i]->getForce(distrt)
                            - offsets[i] * dweights[i];
                  // ffactor += - offsets[i] * dweights[i];
              return ffactor;
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
        RealND cvmu = pot.getColVarMus();
        RealND cvsd = pot.getColVarSds();
        real rc = pot.getCutoff();
        real alp = pot.getAlpha();
        return boost::python::make_tuple(dim, itp, fns, cvrefs,
                                         cvmu, cvsd, alp, rc);
      }
    };

  }
}

#endif
