// ESPP_CLASS
#ifndef _INTERACTION_AKIMA_HPP
#define _INTERACTION_AKIMA_HPP

#include "Interpolation.hpp"



namespace espresso {
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
}//ns espresso




#endif
