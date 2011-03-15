// ESPP_CLASS
#ifndef _INTERACTION_TABULATEDDIHEDRAL_HPP
#define _INTERACTION_TABULATEDDIHEDRAL_HPP

#include "DihedralPotential.hpp"
#include "Interpolation.hpp"

namespace espresso {
    namespace interaction {
     
        class TabulatedDihedral: public DihedralPotentialTemplate <TabulatedDihedral> {
         
            private:
                std::string filename;
                shared_ptr <Interpolation> table;
         
            public:
                static void registerPython();
             
                TabulatedDihedral() {
                    //setCutoff(infinity);
                }
             
                TabulatedDihedral(int itype, const char* filename) {
                    setFilename(itype, filename);
                    std::cout << "using tabulated potential " << filename << "\n";
                }
             
                TabulatedDihedral(int itype, const char* filename, real cutoff) {
                    setFilename(itype, filename);
                    setCutoff(cutoff);
                    std::cout << "using tabulated potential " << filename << "\n";
                }
             
                void setFilename(int itype, const char* _filename);
             
                const char* getFilename() const {
                    return filename.c_str();
                }
             
                real _computeEnergyRaw(real phi) const {
                    if (table)
                        return table->getEnergy(phi);
                    else
                        throw std::runtime_error("Tabulated dihedral potential table not available.");
                }
             
                void _computeForceRaw(Real3D& force1,
                                        Real3D& force2,
                                        Real3D& force3,
                                        Real3D& force4,
                                        const Real3D& dist21,
                                        const Real3D& dist32,
                                        const Real3D& dist43) const {
                    if (table) {
                        // compute phi
                        real dist21_sqr = dist21 * dist21;
                        real dist32_sqr = dist32 * dist32;
                        real dist43_sqr = dist43 * dist43;
                        real dist21_magn = sqrt(dist21_sqr);
                        real dist32_magn = sqrt(dist32_sqr);
                        real dist43_magn = sqrt(dist43_sqr);

                        // cos0
                        real sb1 = 1.0 / dist21_sqr;
                        real sb2 = 1.0 / dist32_sqr;
                        real sb3 = 1.0 / dist43_sqr;
                        real rb1 = sqrt(sb1);
                        real rb3 = sqrt(sb3);
                        real c0 = dist21 * dist43 * rb1 * rb3;


                        // 1st and 2nd angle
                        real ctmp = dist21 * dist32;
                        real r12c1 = 1.0 / (dist21_magn * dist32_magn);
                        real c1mag = ctmp * r12c1;

                        ctmp = (-1.0 * dist32) * dist43;
                        real r12c2 = 1.0 / (dist32_magn * dist43_magn);
                        real c2mag = ctmp * r12c2;


                        //cos and sin of 2 angles and final cos
                        real sin2 = 1.0 - c1mag * c1mag;
                        if (sin2 < 0) sin2 = 0.0;
                        real sc1 = sqrt(sin2);
                        sc1 = 1.0 / sc1;

                        sin2 = 1.0 - c2mag * c2mag;
                        if (sin2 < 0) sin2 = 0.0;
                        real sc2 = sqrt(sin2);
                        sc2 = 1.0 / sc2;

                        real s1 = sc1 * sc1;
                        real s2 = sc2 * sc2;
                        real s12 = sc1 * sc2;
                        real c = (c0 + c1mag * c2mag) * s12;

                        Real3D cc = dist21.cross(dist32);
                        real cmag = sqrt(cc * cc);
                        real dx = cc * dist43 / cmag / dist43_magn;

                        if (c > 1.0) c = 1.0;
                        else if (c < -1.0) c = -1.0;

                        // phi
                        real phi = acos(c);
                        if (dx < 0.0) phi *= -1.0;

                        //phi = 1.0; //testing

                        // read table
                        real a = table->getForce(phi);

                        c = c * a;
                        s12 = s12 * a;

                        real a11 = c * sb1 * s1;
                        real a22 = -sb2 * (2.0 * c0 * s12 - c * (s1 + s2));
                        real a33 = c * sb3 * s2;
                        real a12 = -r12c1 * (c1mag * c * s1 + c2mag * s12);
                        real a13 = -rb1 * rb3 * s12;
                        real a23 = r12c2 * (c2mag * c * s2 + c1mag * s12);

                        Real3D sf2 = a12 * dist21 + a22 * dist32 + a23 * dist43;

                        force1 = a11 * dist21 + a12 * dist32 + a13 * dist43;
                        force2 = (-1.0 * sf2) - force1;
                        force4 = a13 * dist21 + a23 * dist32 + a33 * dist43;
                        force3 = sf2 - force4;
                    }
                    else {
                        throw std::runtime_error("Tabulated dihedral potential table not available.");
                    }
                    
                }
             
                real _computeForceRaw(real phi) const {
                    if (table)
                        return table->getForce(phi);
                    else
                        throw std::runtime_error("Tabulated dihedral potential table not available.");
                }
             
        }; // class
     
    } // ns interaction

} //ns espresso

#endif
