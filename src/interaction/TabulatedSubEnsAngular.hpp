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
#ifndef _INTERACTION_TABULATEDSUBENSANGULAR_HPP
#define _INTERACTION_TABULATEDSUBENSANGULAR_HPP

#include "AngularPotential.hpp"
#include "Interpolation.hpp"
#include "RealND.hpp"
#include "bc/BC.hpp"

namespace espressopp {
    namespace interaction {

        class TabulatedSubEnsAngular: public AngularPotentialTemplate<TabulatedSubEnsAngular> {

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

                TabulatedSubEnsAngular() :
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

                void computeColVarWeights(const Real3D& dist12,
                    const Real3D& dist32, const bc::BC& bc);

                void setColVar(const Real3D& dist12,
                    const Real3D& dist32, const bc::BC& bc);

                real _computeEnergyRaw(real theta) const {
                    real e = 0.;
                    for	(int i=0; i<numInteractions; ++i)
                        e += weights[i] * tables[i]->getEnergy(theta);
                    return e;
                }

                bool _computeForceRaw(Real3D& force12, Real3D& force32,
                                      const Real3D& dist12, const Real3D& dist32) const {
                    real dist12_sqr = dist12 * dist12;
                    real dist32_sqr = dist32 * dist32;
                    real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
                    real cos_theta = dist12 * dist32 / dist1232;
                    real theta = acos(cos_theta);

                    real a = 0.;
                    for	(int i=0; i<numInteractions; ++i)
                        a += weights[i] * tables[i]->getForce(theta);

                    a*=1.0/(sqrt(1.0-cos_theta*cos_theta));

                    real a11 = a * cos_theta / dist12_sqr;
                    real a12 = -a / dist1232;
                    real a22 = a * cos_theta / dist32_sqr;

                    force12 = a11 * dist12 + a12 * dist32;
                    force32 = a22 * dist32 + a12 * dist12;
                    return true;
                }

                real _computeForceRaw(real theta) const {
                    real f = 0.;
                    for	(int i=0; i<numInteractions; ++i)
                        f += weights[i] * tables[i]->getForce(theta);
                    return f;
                }

        }; // class

        // provide pickle support
        struct TabulatedSubEnsAngular_pickle : boost::python::pickle_suite {
            static boost::python::tuple getinitargs(TabulatedSubEnsAngular const& pot) {
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

    } // ns interaction

} //ns espressopp

#endif
