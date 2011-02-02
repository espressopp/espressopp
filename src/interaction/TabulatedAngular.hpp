// ESPP_CLASS
#ifndef _INTERACTION_TABULATEDANGULAR_HPP
#define _INTERACTION_TABULATEDANGULAR_HPP

#include "AngularPotential.hpp"
#include "InterpolationTable.hpp"

namespace espresso {
    namespace interaction {
     
        class TabulatedAngular: public AngularPotentialTemplate <TabulatedAngular> {
         
            private:
                std::string filename;
                shared_ptr <InterpolationTable> table;
         
            public:
                static void registerPython();
             
                TabulatedAngular() {
                    //setCutoff(infinity);
                }
             
                TabulatedAngular(const char* filename) {
                    setFilename(filename);
                }
             
                TabulatedAngular(const char* filename, real cutoff) {
                    setFilename(filename);
                    setCutoff(cutoff);
                }
             
                void setFilename(const char* _filename);
             
                const char* getFilename() const {
                    return filename.c_str();
                }
             
                real _computeEnergyRaw(real theta) const {
                    return table->getEnergy(theta);
                }
             
                void _computeForceRaw(Real3D& force12, Real3D& force32,
                                      const Real3D& dist12, const Real3D& dist32) const {
                    real dist12_sqr = dist12 * dist12;
                    real dist32_sqr = dist32 * dist32;
                    //real dist12_magn = sqrt(dist12_sqr);
                    //real dist32_magn = sqrt(dist32_sqr);
                    real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
                    real cos_theta = dist12 * dist32 / dist1232;
                    //real sin_theta = sqrt(1.0 - cos_theta * cos_theta);
                 
                    real force = table->getForce(acos(cos_theta));
                 
                    real a11 = force * cos_theta / dist12_sqr;
                    //real a12 = -force / (sqrt(dist12_sqr) * sqrt(dist32_sqr));
                    real a12 = -force / dist1232;
                    real a22 = force * cos_theta / dist32_sqr;
                 
                    force12 = a11 * dist12 + a12 * dist32;
                    force32 = a22 * dist32 + a12 * dist12;
                 
                }
             
                real _computeForceRaw(real theta) const {
                    real force = table->getForce(theta);
                    return force;
                }
             
        }; // class
     
    } // ns interaction

} //ns espresso

#endif