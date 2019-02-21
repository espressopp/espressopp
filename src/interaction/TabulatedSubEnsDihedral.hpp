/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research

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
#ifndef _INTERACTION_TABULATEDSUBENSDIHEDRAL_HPP
#define _INTERACTION_TABULATEDSUBENSDIHEDRAL_HPP

#include "DihedralPotential.hpp"
#include "Interpolation.hpp"
#include "RealND.hpp"
#include "bc/BC.hpp"

namespace espressopp {
    namespace interaction {

        class TabulatedSubEnsDihedral: public DihedralPotentialTemplate <TabulatedSubEnsDihedral> {

            private:
                int numInteractions;
                std::vector<std::string> filenames;
                std::vector<shared_ptr <Interpolation>> tables;
                int interpolationType;
                // Reference values of the collective variable centers
                RealNDs colVarRef;
                // Weights of each table
                RealND weights;
                // Target probability of each table
                RealND targetProb;
                // Running sum of each weight and number of counts
                RealND weightSum;
                int weightCounts;
                // Renormalize collective variables: std
                RealND colVarSd;
                // characteristic decay length of the interpolation
                real alpha;
                // Size of CV partners
                int colVarBondListSize;
                int colVarAngleListSize;
                int colVarDihedListSize;

            public:
                static void registerPython();

                TabulatedSubEnsDihedral() :
                    numInteractions(0) {
                    setCutoff(infinity);
                    weights.setDimension(0);
                    weightSum.setDimension(0);
                    targetProb.setDimension(0);
                    weightCounts = 0;
                    colVarSd.setDimension(3);
                    colVarRef.setDimension(0);
                    alpha = 1.;
                    colVarBondListSize = 0;
                    colVarAngleListSize = 0;
                    colVarDihedListSize = 0;
                }

                void addInteraction(int itype, boost::python::str fname,
                                    const RealND& _cvref);

                void setDimension(int _dim) {
                  numInteractions = _dim;
                  colVarRef.setDimension( numInteractions );
                  tables.resize( numInteractions );
                  filenames.resize( numInteractions );
                  weights.setDimension( numInteractions );
                  weightSum.setDimension( numInteractions );
                  targetProb.setDimension( numInteractions );
                }

                int getDimension() const { return numInteractions; }

                /** Setter for the interpolation type */
                void setInterpolationType(int itype) { interpolationType = itype; }

                /** Getter for the interpolation type */
                int getInterpolationType() const { return interpolationType; }

                RealND getTargetProb() const { return targetProb; }

                void setTargetProb(int index, real _r) {
                    return targetProb.setItem(index, _r);
                }

                RealND getColVarSds() const { return colVarSd; }

                void setColVarSd(int index, real _r) {
                    return colVarSd.setItem(index, _r);
                }

                RealND getColVarRef(int i) const { return colVarRef[i]; }

                RealNDs getColVarRefs() const { return colVarRef; }

                void setColVarRef(const RealNDs& cvRefs);

                void setColVarRefs(const RealNDs& c) { colVarRef = c; }

                boost::python::list getFilenames() const {
                    return boost::python::list(filenames); }

                boost::python::list getFilename(int index) const {
                    return boost::python::list(filenames[index]);
                }

                void setFilename(int index, boost::python::list _f) {
                    filenames[index] = boost::python::extract<std::string>(_f);
                }

                void setFilenames(int dim, int itype, boost::python::list _filenames);

                RealND getWeights() const { return weights; }

                void setWeights(const RealND& r) { weights = r; }

                real getWeight(int index) const {
                    return weights.getItem(index);
                }

                void setWeight(int index, real _w) {
                    return weights.setItem(index, _w);
                }

                real getAlpha() const { return alpha; }

                void setAlpha(real _r) { alpha = _r; }

                void computeColVarWeights(const Real3D& dist21, const Real3D& dist32,
                                const Real3D& dist43, const bc::BC& bc);

                void setColVar(const Real3D& dist21, const Real3D& dist32,
                                const Real3D& dist43, const bc::BC& bc);

                real _computeEnergyRaw(real phi) const {
                    real e = 0.;
                    for	(int i=0; i<numInteractions; ++i)
                        e += weights[i] * tables[i]->getEnergy(phi);
                    return e;
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
                    real a = 0.;
                    for	(int i=0; i<numInteractions; ++i)
                        a += weights[i] * tables[i]->getForce(_phi);

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

                real _computeForceRaw(real phi) const {
                    real f = 0.;
                    for	(int i=0; i<numInteractions; ++i)
                        f += weights[i] * tables[i]->getForce(phi);
                    return f;
                }

        }; // class
        // provide pickle support
        struct TabulatedSubEnsDihedral_pickle : boost::python::pickle_suite
        {
          static
          boost::python::tuple
          getinitargs(TabulatedSubEnsDihedral const& pot)
          {
              int itp = pot.getInterpolationType();
              boost::python::list fns;
              RealNDs cvrefs = pot.getColVarRefs();
              int dim = pot.getDimension();
              fns = pot.getFilenames();
              RealND cvsd = pot.getColVarSds();
              real rc = pot.getCutoff();
              real alp = pot.getAlpha();
              return boost::python::make_tuple(dim, itp, fns, cvrefs,
                                                cvsd, alp, rc);
          }
        };

    } //ns espressopp

}
#endif
