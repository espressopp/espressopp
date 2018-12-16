/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research
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
#ifndef _INTERACTION_TABULATEDDIHEDRAL_HPP
#define _INTERACTION_TABULATEDDIHEDRAL_HPP

#include "DihedralPotential.hpp"
#include "Interpolation.hpp"

namespace espressopp {
    namespace interaction {

        class TabulatedDihedral: public DihedralPotentialTemplate <TabulatedDihedral> {

            private:
                std::string filename;
                shared_ptr <Interpolation> table;
                int interpolationType;

            public:
                static void registerPython();

                TabulatedDihedral() {
                    //setCutoff(infinity);
                }

                TabulatedDihedral(int itype, const char* filename) {
                    setFilename(itype, filename);
                    setInterpolationType(itype);
                }

                TabulatedDihedral(int itype, const char* filename, real cutoff) {
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

                real _computeEnergyRaw(real phi) const {
                    if (table)
                        return table->getEnergy(phi);
                    else
                        throw std::runtime_error("Tabulated dihedral potential table not available.");
                }

                // Kroneker delta function
                real d(int i, int j) const {
                  if(i==j)
                    return 1.0;
                  else
                    return 0.0;
                }

                /** Compute dot product of two vectors */
                Real3D prod(Real3D a, Real3D b) const {
                  Real3D res(0.0, 0.0, 0.0);
                  for(int i=0; i<3; i++){
                    for(int j=0; j<3; j++){
                      res[i]+=(1.0-d(i, j))*a[j]*b[j];
                    }
                  }
                  return res;
                }

                void _computeForceRaw(Real3D& force1,
                                        Real3D& force2,
                                        Real3D& force3,
                                        Real3D& force4,
                                        const Real3D& r21,
                                        const Real3D& r32,
                                        const Real3D& r43) const {
                    if (table) {
                      Real3D retF[4];

                      Real3D rijjk = r21.cross(r32); // [r21 x r32]
                      Real3D rjkkn = r32.cross(r43); // [r32 x r43]

                      real rijjk_sqr = rijjk.sqr();
                      real rjkkn_sqr = rjkkn.sqr();

                      real rijjk_abs = sqrt(rijjk_sqr);
                      real rjkkn_abs = sqrt(rjkkn_sqr);

                      real inv_rijjk = 1.0 / rijjk_abs;
                      real inv_rjkkn = 1.0 / rjkkn_abs;

                      // cosine between planes
                      real cos_phi = (rijjk * rjkkn) * (inv_rijjk * inv_rjkkn);
                      real _phi = acos(cos_phi);
                      if (cos_phi > 1.0) {
                        cos_phi = 1.0;
                        _phi = 1e-10; //not 0.0, because 1.0/sin(_phi) would cause a singularity
                      } else if (cos_phi < -1.0) {
                        cos_phi = -1.0;
                        _phi = M_PI-1e-10;
                      }

                      //get sign of phi
                      //positive if (rij x rjk) x (rjk x rkn) is in the same direction as rjk, negative otherwise (see DLPOLY manual)
                      Real3D rcross = rijjk.cross(rjkkn); //(rij x rjk) x (rjk x rkn)
                      real signcheck = rcross * r32;
                      if (signcheck < 0.0) _phi *= -1.0;

                      // read table
                      real a = table->getForce(_phi);

                      real coef1 = -(1.0/sin(_phi)) * a;

                      real A1 = inv_rijjk * inv_rjkkn;
                      real A2 = inv_rijjk * inv_rijjk;
                      real A3 = inv_rjkkn * inv_rjkkn;

                      Real3D p3232 ( prod(r32,r32) );
                      Real3D p3243 ( prod(r32,r43) );
                      Real3D p2132 ( prod(r21,r32) );
                      Real3D p2143 ( prod(r21,r43) );
                      Real3D p2121 ( prod(r21,r21) );
                      Real3D p4343 ( prod(r43,r43) );

                      // we have 4 particles 1,2,3,4
                      for(int l=0; l<4; l++){
                        Real3D B1, B2, B3;

                        for(int i=0; i<3; i++){
                          B1[i]= r21[i] * ( p3232[i] * (d(l,2)-d(l,3)) + p3243[i] * (d(l,2)-d(l,1)) ) +
                                 r32[i] * ( p2132[i] * (d(l,3)-d(l,2)) + p3243[i] * (d(l,1)-d(l,0)) ) +
                                 r43[i] * ( p2132[i] * (d(l,2)-d(l,1)) + p3232[i] * (d(l,0)-d(l,1)) ) +
                                 2.0 * r32[i] * p2143[i] * (d(l,1)-d(l,2));

                          B2[i]= 2.0 * r21[i] * ( p3232[i] * (d(l,1)-d(l,0)) + p2132[i] * (d(l,1)-d(l,2)) ) +
                                 2.0 * r32[i] * ( p2121[i] * (d(l,2)-d(l,1)) + p2132[i] * (d(l,0)-d(l,1)));

                          B3[i]= 2.0 * r43[i] * ( p3232[i] * (d(l,3)-d(l,2)) + p3243[i] * (d(l,1)-d(l,2)) ) +
                                 2.0 * r32[i] * ( p4343[i] * (d(l,2)-d(l,1)) + p3243[i] * (d(l,2)-d(l,3)));
                        }

                        retF[l] = coef1 * ( A1*B1 - 0.5*cos_phi*(A2*B2 + A3*B3) );
                      }

                      force1 = retF[0];
                      force2 = retF[1];
                      force3 = retF[2];
                      force4 = retF[3];

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
    // provide pickle support
    struct TabulatedDihedral_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(TabulatedDihedral const& pot)
      {   int itp        = pot.getInterpolationType();
          std::string fn = pot.getFilename();
          real rc        = pot.getCutoff();
          return boost::python::make_tuple(itp, fn,rc);
      }
    };

    } // ns interaction

} //ns espressopp

#endif
