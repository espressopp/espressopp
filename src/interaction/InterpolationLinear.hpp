// ESPP_CLASS
#ifndef _INTERACTION_LINEAR_HPP
#define _INTERACTION_LINEAR_HPP

#include "Interpolation.hpp"



namespace espresso {
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
}//ns espresso




#endif