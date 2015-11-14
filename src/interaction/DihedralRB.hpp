/*
  Copyright (C) 2015
      Jakub Krajniak (jkrajniak at gmail.com)

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
#ifndef _INTERACTION_DihedralRB_HPP
#define _INTERACTION_DihedralRB_HPP

#include <cmath>
#include "DihedralPotential.hpp"

namespace espressopp {
namespace interaction {

/** This class provides methods to compute forces and energies of
    the Ryckaert-Bellemans  dihedral potential.*/
class DihedralRB : public DihedralPotentialTemplate< DihedralRB > {
 private:
  real K0;
  real K1;
  real K2;
  real K3;
  real K4;
  real K5;

 public:
  static void registerPython();

  DihedralRB() : K0(0.0), K1(0.0), K2(0.0), K3(0.0), K4(0.0), K5(0.0) { }
  DihedralRB(real _K0, real _K1, real _K2, real _K3, real _K4, real _K5)
      : K0(_K0), K1(_K1), K2(_K2), K3(_K3), K4(_K4), K5(_K5) { }

  real _computeEnergyRaw(real _phi) const {
    // _phi should be in radians
    real cos_phi = cos(_phi);
    real energy = K0 +
        K1 * cos_phi +
        K2 * cos_phi*cos_phi +
        K3 * cos_phi*cos_phi*cos_phi +
        K4 * cos_phi*cos_phi*cos_phi*cos_phi +
        K5*cos_phi*cos_phi*cos_phi*cos_phi*cos_phi;
    return energy;
  }

  // Kroneker delta function
  real d(int i, int j) const {
    if (i == j)
      return 1.0;
    else
      return 0.0;
  }

  Real3D prod(Real3D a, Real3D b) const {
    Real3D res(0.0, 0.0, 0.0);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        res[i] += (1.0-d(i, j))*a[j]*b[j];
    return res;
  }

  void _computeForceRaw(Real3D& force1,
                        Real3D& force2,
                        Real3D& force3,
                        Real3D& force4,
                        const Real3D& r21,
                        const Real3D& r32,
                        const Real3D& r43) const {
    Real3D retF[4];

    Real3D rijjk = r21.cross(r32);  // [r21 x r32]
    Real3D rjkkn = r32.cross(r43);  // [r32 x r43]

    real rijjk_sqr = rijjk.sqr();
    real rjkkn_sqr = rjkkn.sqr();

    real rijjk_abs = sqrt(rijjk_sqr);
    real rjkkn_abs = sqrt(rjkkn_sqr);

    // If vectors are collinear then return, it is not
    // possible to create a plane and define an angle.
    if (rijjk_abs < 0.00000001 || rjkkn_abs < 0.00000001)
      return;

    real inv_rijjk = 1.0 / rijjk_abs;
    real inv_rjkkn = 1.0 / rjkkn_abs;

    // cosine between planes
    real cos_phi = (rijjk * rjkkn) * (inv_rijjk * inv_rjkkn);
    if (cos_phi > 1.0) cos_phi = 1.0;
    else if (cos_phi < -1.0) cos_phi = -1.0;

    /** Derivative of the potential. */
    real coef1 = -K1
      - 2.0*K2*cos_phi
      - 3.0*K3*cos_phi*cos_phi
      - 4.0*K4*cos_phi*cos_phi*cos_phi
      - 5.0*K5*cos_phi*cos_phi*cos_phi*cos_phi;

    /** Calculates force in Cartesian coordinates.
     * base on: http://www.ccp5.ac.uk/DL_POLY_CLASSIC/MANUALS/USRMAN.pdf
     */
    real A1 = inv_rijjk * inv_rjkkn;
    real A2 = inv_rijjk * inv_rijjk;
    real A3 = inv_rjkkn * inv_rjkkn;

    Real3D p3232(prod(r32, r32));
    Real3D p3243(prod(r32, r43));
    Real3D p2132(prod(r21, r32));
    Real3D p2143(prod(r21, r43));
    Real3D p2121(prod(r21, r21));
    Real3D p4343(prod(r43, r43));

    // we have 4 particles 1,2,3,4
    for (int l = 0; l < 4; l++) {
      Real3D B1, B2, B3;

      for (int i = 0; i < 3; i++) {
        B1[i] = r21[i] * (p3232[i] * (d(l, 2) - d(l, 3)) + p3243[i] * (d(l, 2) - d(l, 1))) +
                r32[i] * (p2132[i] * (d(l, 3) - d(l, 2)) + p3243[i] * (d(l, 1) - d(l, 0))) +
                r43[i] * (p2132[i] * (d(l, 2) - d(l, 1)) + p3232[i] * (d(l, 0) - d(l, 1))) +
                2.0 * r32[i] * p2143[i] * (d(l, 1) - d(l, 2));

        B2[i] = 2.0 * r21[i] * (p3232[i] * (d(l, 1) - d(l, 0)) + p2132[i] * (d(l, 1) - d(l, 2))) +
                2.0 * r32[i] * (p2121[i] * (d(l, 2) - d(l, 1)) + p2132[i] * (d(l, 0) - d(l, 1)));

        B3[i] = 2.0 * r43[i] * (p3232[i] * (d(l, 3) - d(l, 2)) + p3243[i] * (d(l, 1) - d(l, 2))) +
                2.0 * r32[i] * (p4343[i] * (d(l, 2) - d(l, 1)) + p3243[i] * (d(l, 2) - d(l, 3)));
      }

      retF[l] = coef1 * (A1*B1 - 0.5*cos_phi*(A2*B2 + A3*B3));
    }

    force1 = retF[0];
    force2 = retF[1];
    force3 = retF[2];
    force4 = retF[3];
  }

  // used for generating tabulated potential file
  real _computeForceRaw(real phi) const {
    real cos_phi = cos(phi);
    real coef1 = -K1
      - 2.0*K2*cos_phi
      - 3.0*K3*cos_phi*cos_phi
      - 4.0*K4*cos_phi*cos_phi*cos_phi
      - 5.0*K5*cos_phi*cos_phi*cos_phi*cos_phi;
    return -1.0*coef1;
  }
};

}  // end namespace interaction
}  // end namespace espressopp

#endif
