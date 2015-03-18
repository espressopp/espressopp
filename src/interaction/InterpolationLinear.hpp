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
#ifndef _INTERACTION_LINEAR_HPP
#define _INTERACTION_LINEAR_HPP

#include "Interpolation.hpp"



namespace espressopp {
    namespace interaction {
        class InterpolationLinear: public InterpolationTemplate <InterpolationLinear> {
            public:
                InterpolationLinear();
                ~InterpolationLinear();
                void readRaw(mpi::communicator comm, const char* file);
                real getEnergyRaw(real r) const;
                real getForceRaw(real r) const;
            
            protected:
                static LOG4ESPP_DECL_LOGGER(theLogger);
            
            private:
                // make copy constructor private because copy is not allowed
                InterpolationLinear(const InterpolationLinear&) {}
                
                // Reading values from file, control processor only; returns
                //  number of valid entries, error if value is less than 2.
                // If dummy is true, values will not be stored in arrays r, e, f
                int readFile(const char* file, bool dummy);
                
                // Spline read-in values
                void spline(const real* x, const real* y, int N, real* a, real* b);
                
                // Spline interpolation
                real splineInterpolation(real r, const real* a, const real* b) const;
             
                int N;  // number of read values
             
                real inner;
                real outer;
                real delta;
                real invdelta;
             
                real *radius, *energy, *force;
                
                real *ae, *be;
                real *af, *bf;
            
        };//class InterpolationLinear
        
        
        inline real InterpolationLinear::getEnergyRaw(real r) const {
            return splineInterpolation(r, ae, be);
        }
        
        inline real InterpolationLinear::getForceRaw(real r) const {
            return splineInterpolation(r, af, bf);
        }
        
        inline real InterpolationLinear::splineInterpolation(real r,
                            const real* a, const real* b) const {
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
            
            return a[index]*r + b[index];
        }
        
        
    }//ns interaction
}//ns espressopp




#endif