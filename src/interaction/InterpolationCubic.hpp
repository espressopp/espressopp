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
#ifndef _INTERACTION_CUBIC_HPP
#define _INTERACTION_CUBIC_HPP

#include "Interpolation.hpp"



namespace espressopp {
    namespace interaction {
        class InterpolationCubic: public InterpolationTemplate <InterpolationCubic> {
            public:
                InterpolationCubic();
                ~InterpolationCubic();
                void readRaw(mpi::communicator comm, const char* file);
                real getEnergyRaw(real r) const;
                real getForceRaw(real r) const;
            
            protected:
                static LOG4ESPP_DECL_LOGGER(theLogger);
            
            private:
                // make copy constructor private because copy is not allowed
                InterpolationCubic(const InterpolationCubic&) {}
                
                // Reading values from file, control processor only; returns
                //  number of valid entries, error if value is less than 2.
                // If dummy is true, values will not be stored in arrays r, e, f
                int readFile(const char* file, bool dummy);
                
                // Spline read-in values
                void spline(const real* x, const real* y, int n,
                            real yp1, real ypn, real* y2);
                
                // Spline interpolation
                real splineInterpolation(real r, const real* fn, const real* fn2) const;
             
                int N;  // number of read values
             
                real inner;
                real outer;
                real delta;
                real invdelta;
                real deltasq6;
             
                bool allocated;
             
                real *radius;
                real *energy;
                real *force;
             
                real *energy2;  // used for spline interpolation
                real *force2;   // used for spline interpolation
            
        };//class InterpolationCubic
        
        
        inline real InterpolationCubic::getEnergyRaw(real r) const {
            return splineInterpolation(r, energy, energy2);
        }
        
        inline real InterpolationCubic::getForceRaw(real r) const {
            return splineInterpolation(r, force, force2);
        }
        
        inline real InterpolationCubic::splineInterpolation(real r, const real* fn, 
                                                         const real* fn2) const {
            int index = static_cast<int>((r - inner) * invdelta);
                
            if (index < 0) {
                LOG4ESPP_ERROR(theLogger, "distance " << r << " out of range "
                                << inner << " - " << inner + (N - 1) * delta);
                index = 0;
            }
            else if (index >= N) {
                LOG4ESPP_ERROR(theLogger, "distance " << r << " out of range "
                                << inner << " - " << inner + (N - 1) * delta);
                index = N-1;
            }
                
            real b = (r - radius[index]) * invdelta;
            real a = 1.0 - b;
            real f = a * fn[index] +
                    b * fn[index+1] +
                    ((a*a*a-a)*fn2[index] +
                        (b*b*b-b)*fn2[index+1]) *
                    deltasq6;
                
            // printf("interpoloate %f, a = %f, b = %f, fn = %f - %f, fn2 = %f - %f\n", 
            //        r, a, b, fn[index], fn[index+1], fn2[index], fn2[index+1]);
            return f;
        }
        
        
    }//ns interaction
}//ns espressopp




#endif