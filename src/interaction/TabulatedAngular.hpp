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
#ifndef _INTERACTION_TABULATEDANGULAR_HPP
#define _INTERACTION_TABULATEDANGULAR_HPP

#include "AngularPotential.hpp"
#include "Interpolation.hpp"

namespace espressopp {
    namespace interaction {
     
        class TabulatedAngular: public AngularPotentialTemplate <TabulatedAngular> {
         
            private:
                std::string filename;
                shared_ptr <Interpolation> table;
                int interpolationType;
         
            public:
                static void registerPython();
             
                TabulatedAngular() {
                    //setCutoff(infinity);
                    //std::cout << "using default tabulated potential ...\n";
                }
             
                TabulatedAngular(int itype, const char* filename) {
                    setFilename(itype, filename);
                    setInterpolationType(itype);
                }
             
                TabulatedAngular(int itype, const char* filename, real cutoff) {
                    setFilename(itype, filename);
                    setInterpolationType(itype);
                    setCutoff(cutoff);
                    std::cout << "using tabulated potential " << filename << "\n";
                }
                /** Setter for the interpolation type */
                void setInterpolationType(int itype) { interpolationType = itype; }

                /** Getter for the interpolation type */
                int getInterpolationType() const { return interpolationType; }

                void setFilename(int itype, const char* _filename);
             
                const char* getFilename() const {
                    return filename.c_str();
                }
             
                real _computeEnergyRaw(real theta) const {
                    if (table) {
                      return table->getEnergy(theta);
                    } else {
                      LOG4ESPP_DEBUG(theLogger, "Tabulate angular potential table not available.");
                      return 0.0;
                    }
                }
             
                bool _computeForceRaw(Real3D& force12, Real3D& force32,
                                      const Real3D& dist12, const Real3D& dist32) const {
                    if (table) {
                        real dist12_sqr = dist12 * dist12;
                        real dist32_sqr = dist32 * dist32;
                        real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
                        real cos_theta = dist12 * dist32 / dist1232;

                        real a = table->getForce(acos(cos_theta));
                        
                        a*=1.0/(sqrt(1.0-cos_theta*cos_theta));
                        
                        real a11 = a * cos_theta / dist12_sqr;
                        real a12 = -a / dist1232;
                        real a22 = a * cos_theta / dist32_sqr;

                        force12 = a11 * dist12 + a12 * dist32;
                        force32 = a22 * dist32 + a12 * dist12;
                    }
                    else {
                        LOG4ESPP_DEBUG(theLogger, "Tabulate angular potential table not available.");
                        return false;
                    }
                    return true;
                }
             
                real _computeForceRaw(real theta) const {
                    return table->getForce(theta);
                }
             
        }; // class

        // provide pickle support
        struct TabulatedAngular_pickle : boost::python::pickle_suite {
          static boost::python::tuple getinitargs(TabulatedAngular const& pot) {
            int itp = pot.getInterpolationType();
            std::string fn = pot.getFilename();
            real rc = pot.getCutoff();
            return boost::python::make_tuple(itp, fn,rc);
          }
        };
     
    } // ns interaction

} //ns espressopp

#endif
