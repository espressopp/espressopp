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
#ifndef _INTERACTION_AKIMA_HPP
#define _INTERACTION_AKIMA_HPP

#include "Interpolation.hpp"



namespace espressopp {
    namespace interaction {
        class InterpolationAkima: public InterpolationTemplate <InterpolationAkima> {
            public:
                InterpolationAkima();
                ~InterpolationAkima();
                void readRaw(mpi::communicator comm, const char* file);
                real getEnergyRaw(real r) const;
                real getForceRaw(real r) const;
            
            protected:
                static LOG4ESPP_DECL_LOGGER(theLogger);
            
            private:
                // make copy constructor private because copy is not allowed
                InterpolationAkima(const InterpolationAkima&) {}
                
                // Reading values from file, control processor only; returns
                //  number of valid entries, error if value is less than 2.
                // If dummy is true, values will not be stored in arrays r, e, f
                int readFile(const char* file, bool dummy);
                
                // Spline read-in values
                void spline(const real* x, const real* y, int N,
                                    real* p0, real* p1, real* p2, real* p3);
                
                // Spline interpolation
                real splineInterpolation(real r,
                            const real* p0, const real* p1, const real* p2, const real* p3) const;
                real getSlope(real m1, real m2, real m3, real m4);
             
                int N;  // number of read values
             
                real inner;
                real outer;
                real delta;
                real invdelta;
             
                real *radius, *energy, *force;
                
                real *p0e, *p1e, *p2e, *p3e;
                real *p0f, *p1f, *p2f, *p3f;
            
        };//class InterpolationAkima
        
        
        inline real InterpolationAkima::getEnergyRaw(real r) const {
            return splineInterpolation(r, p0e, p1e, p2e, p3e);
        }
        
        inline real InterpolationAkima::getForceRaw(real r) const {
            return splineInterpolation(r, p0f, p1f, p2f, p3f);
        }
        
        inline real InterpolationAkima::splineInterpolation(real r,
                            const real* p0, const real* p1, const real* p2, const real* p3) const {
            int index;
            index = static_cast<int>((r - inner) * invdelta);
            if (index < 0) {
                LOG4ESPP_ERROR(theLogger, "distance " << r << " out of range "
                                << inner << " - " << inner + (N - 1) * delta);
                index = 0;
            }
            if (index >= N) {
                LOG4ESPP_ERROR(theLogger, "distance " << r << " out of range "
                                << inner << " - " << inner + (N - 1) * delta);
                index = N-1;
            }
            
            real z = r - radius[index];
            int zz2 = z*z;
            return p0[index] +
                   p1[index] * z +
                   p2[index] * zz2 +
                   p3[index] * zz2 * z;
        }
        
        inline real InterpolationAkima::getSlope(real m1, real m2, real m3, real m4) {
            if ((m1 == m2) && (m3 == m4)) {
                return (m2 + m3) / 2.0;
            } else {
                return (fabs(m4 - m3) * m2 + fabs(m2 - m1) * m3) /
                        (fabs(m4 - m3) + fabs(m2 - m1));
            }
        }
    }//ns interaction
}//ns espressopp




#endif
