// ESPP_CLASS
#ifndef _INTERACTION_TABULATEDANGULAR_HPP
#define _INTERACTION_TABULATEDANGULAR_HPP

#include "AngularPotential.hpp"
#include "Interpolation.hpp"

namespace espresso {
    namespace interaction {
     
        class TabulatedAngular: public AngularPotentialTemplate <TabulatedAngular> {
         
            private:
                std::string filename;
                shared_ptr <Interpolation> table;
         
            public:
                static void registerPython();
             
                TabulatedAngular() {
                    //setCutoff(infinity);
                }
             
                TabulatedAngular(int itype, const char* filename) {
                    setFilename(itype, filename);
                }
             
                TabulatedAngular(int itype, const char* filename, real cutoff) {
                    setFilename(itype, filename);
                    setCutoff(cutoff);
                }
             
                void setFilename(int itype, const char* _filename);
             
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
                    real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
                    real cos_theta = dist12 * dist32 / dist1232;
                 
                    real a = table->getForce(acos(cos_theta));
                 
                    real a11 = a * cos_theta / dist12_sqr;
                    real a12 = -a / dist1232;
                    real a22 = a * cos_theta / dist32_sqr;
                 
                    force12 = a11 * dist12 + a12 * dist32;
                    force32 = a22 * dist32 + a12 * dist12;
                 
                }
             
                real _computeForceRaw(real theta) const {
                    return table->getForce(theta);
                }
             
        }; // class
     
    } // ns interaction

} //ns espresso

#endif