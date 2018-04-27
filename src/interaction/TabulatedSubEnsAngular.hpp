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

// ESPP_CLASS
#ifndef _INTERACTION_TABULATEDSUBENSANGULAR_HPP
#define _INTERACTION_TABULATEDSUBENSANGULAR_HPP

#include "AngularPotential.hpp"
#include "Interpolation.hpp"

namespace espressopp {
    namespace interaction {

        class TabulatedSubEnsAngular: public AngularPotentialTemplate <TabulatedSubEnsAngular> {

            private:
                int numInteractions;
                std::vector<std::string> filenames;
                std::vector<shared_ptr <Interpolation>> tables;
                int interpolationType;
                // Reference values of the collective variable centers
                std::vector<cv_angle> colVarRef;
                // Weights of each table
                std::vector<real> weights;
                // Scaling factor for the weight
                real alpha = 1.0;

            public:
                static void registerPython();

                TabulatedSubEnsAngular() {
                    //setCutoff(infinity);
                    //std::cout << "using default tabulated potential ...\n";
                }

                TabulatedSubEnsAngular(int dim, int itype, std::vector<std::string> filenames) {
                    setDimension(dim);
                    setFilenames(dim, itype, filenames);
                    setInterpolationType(itype);
                    colVarRef.resize(dim);
                    weights.resize(dim);
                }

                TabulatedSubEnsAngular(int dim, int itype, std::vector<std::string> filenames, real cutoff) {
                    setDimension(dim);
                    setFilenames(dim, itype, filenames);
                    setInterpolationType(itype);
                    colVarRef.resize(dim);
                    weights.resize(dim);
                    setCutoff(cutoff);
                    std::cout << "using tabulated potentials \n";
                    for (int i=0; i<dim; ++i)
                        std::cout << "  " << filenames[i] << "\n";
                }
                void setDimension(int _dim) {
                  numInteractions = _dim;
                  tables.resize( numInteractions );
                }

                int getDimension() const { return numInteractions; }

                /** Setter for the interpolation type */
                void setInterpolationType(int itype) { interpolationType = itype; }

                /** Getter for the interpolation type */
                int getInterpolationType() const { return interpolationType; }

                real getWeightScalingFactor() const { return alpha; }

                void setWeightScalingFactor(real factor) { alpha = factor; }

                cv_angle getColVarRef(int i) const { return colVarRef[i]; }

                void setColVarRef(std::vector< cv_angle > cvRefs);

                double distColVars(cv_angle cv1, cv_angle cv2);

                std::vector<std::string> getFilenames() const { return filenames; }

                void setFilenames(int dim, int itype, std::vector<std::string> _filenames);

                std::vector<real> getWeights() const { return weights; };

                void computeColVarWeights(const Real3D& dist12, const Real3D& dist32);

                real _computeEnergyRaw(real theta) const {
                    real e = 0.;
                    for (int i=0; i<numInteractions; ++i) {
                        if (tables[i]) {
                            e += weights[i] * tables[i]->getEnergy(theta);
                        } else {
                            LOG4ESPP_DEBUG(theLogger,
                                "Tabulated angular potential table not available.");
                            return 0.0;
                        }
                    }
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
                    for (int i=0; i<numInteractions; ++i) {
                        if (tables[i])
                            a += weights[i] * tables[i]->getForce(acos(cos_theta));
                        else {
                            LOG4ESPP_DEBUG(theLogger,
                                "Tabulated angular potential table not available.");
                            return false;
                        }
                    }

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
                    for (int i=0; i<numInteractions; ++i)
                        f += weights[i] * tables[i]->getForce(theta);
                    return f;
                }

        }; // class

        // provide pickle support
        struct TabulatedSubEnsAngular_pickle : boost::python::pickle_suite {
            static boost::python::tuple getinitargs(TabulatedSubEnsAngular const& pot) {
                int itp = pot.getInterpolationType();
                std::vector<std::string> fns;
                std::vector<cv_angle> cvrefs;
                int dim = pot.getDimension();
                fns = pot.getFilenames();
                cvrefs.resize(dim);
                for (int i=0; i<dim; ++i) {
                  cvrefs[i] = pot.getColVarRef(i);
                }
                real alp = pot.getWeightScalingFactor();
                real rc = pot.getCutoff();
                return boost::python::make_tuple(dim, itp, fns, cvrefs, alp, rc);
            }
        };

    } // ns interaction

} //ns espressopp

#endif
