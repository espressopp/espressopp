/*
  Copyright (C) 2012,2013,2016
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
#ifndef _INTERACTION_REACTIONFIELDGENERALIZEDTI_HPP
#define _INTERACTION_REACTIONFIELDGENERALIZEDTI_HPP

#include "Potential.hpp"
#include <algorithm>
#include <set>

using namespace std;

namespace espressopp {
    namespace interaction {
        /*
         * This class provides methods to compute forces and energies of
         * the generalized reaction field
         * in a Thermodynamic Integration simulation
         */
        class ReactionFieldGeneralizedTI: public
        PotentialTemplate<ReactionFieldGeneralizedTI> {
            private:
                real kappa;
                real epsilon1, epsilon2;
                real rc, rc3; // cutoff, cutoff^3
                real rc2;
                real B1, B0, B1_half;
                real prefactor;
                bool annihilate; //if true, atoms in pidsTI are annihilated, otherwise they are decoupled
                real crf; /* const to make potential zero at cutoff, corresponding
                 to the shift, which we cannot use because it depends on qq*/
                real lambdaTI; //not to be confused with the lambda used in AdResS simulations
                real complLambdaTI; //1-lambdaTI
                std::set<longint> pidsTI; //PIDs of particles whose charge is zero in TI state B

                void initialize() {
                    real krc = kappa*rc;
                    real tmp1 = (epsilon1 - 4.0*epsilon2) *
                            (1.0 + krc) - 2.0*epsilon2 * krc*krc;
                    real tmp2 = (epsilon1 + 2.0*epsilon2) *
                            (1.0 + krc) + epsilon2 * krc*krc;
                    rc3 = pow(rc,3);
                    rc2 = pow(rc,2);
                    B1 = tmp1/tmp2;
                    B1 = (1.0+B1) / rc3;
                    B1_half = B1/2.0;
                    crf = 3*epsilon2/(rc*(2*epsilon2+1.0));
                    prefactor/=epsilon1;
                    complLambdaTI = 1.0 - lambdaTI;
                }

            public:
                static void registerPython();

                ReactionFieldGeneralizedTI()
                : prefactor(0.0), kappa(0.0),
                 epsilon1(1.0), epsilon2(80.0),
                 rc(1.0), lambdaTI(0.0), annihilate(1) {
                    setShift(0.0);
                    autoShift = false;
                    setCutoff(infinity);
                    initialize();
                }

                ReactionFieldGeneralizedTI(   
                        real _prefactor,
                        real _kappa,
                        real _eps1,
                        real _eps2,
                        real _cutoff,
                        real _lambdaTI,
                        bool _annihilate)
                : prefactor(_prefactor),kappa(_kappa),
                  epsilon1(_eps1), epsilon2(_eps2),
                  rc(_cutoff), lambdaTI(_lambdaTI), 
                  annihilate(_annihilate) {
                    setShift(0.0);
                    autoShift = false;
                    /*Note: AutoShift cannot be used here since the shift
                     *  has to depend on the product of charges */
                    setCutoff(_cutoff);
                    initialize();
                }

                real getKappa() const { return kappa; }
                real getEpsilon1() const { return epsilon1; }
                real getEpsilon2() const { return epsilon2; }
                real getLambdaTI() const { return lambdaTI; }
                bool getAnnihilate() const {return annihilate; }

                void setCutoff(real _cutoff) {
                    cutoff = _cutoff;
                    cutoffSqr = cutoff*cutoff;
                }

                // Setter and getter
                void setPrefactor(real _prefactor) {
                    prefactor = _prefactor;
                }
                real getPrefactor() const { return prefactor; }

                void addPid(longint pid) {
                  pidsTI.insert(pid);
                }

                bool checkTIpair(size_t pid1, size_t pid2) const {
                  if (annihilate) {
                      //if at least one of pid1 and pid2 in pidsTI (OR)
                      if ((pidsTI.find(pid1) != pidsTI.end()) or (pidsTI.find(pid2) != pidsTI.end())) {
                        return true;
                      } else {
                        return false;
                      }
                  } else { //decouple
                      //if pid1 or pid2 in pidsTI, but not both (XOR)
                      if ((pidsTI.find(pid1) != pidsTI.end()) != (pidsTI.find(pid2) != pidsTI.end())) {
                        return true;
                      } else {
                        return false;
                      }
                  }
                }

                real _computeEnergy(const Particle& p1, const Particle& p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real distSqr = dist.sqr();
                    if (distSqr>rc2) return 0.0;
                    real qq = p1.q()*p2.q();
                    if (checkTIpair(p1.id(),p2.id())) {
                      qq *= complLambdaTI; //(1-lambda)*qAi*qAj
                    }
                    /* Note: this implementation counts minus integral of the force as energy
                         which is not the same as the full electrostatic energy
                         See the original paper Tironi et al, J. Chem. Phys. 102 , 5451 (1995) 
                         for details*/

                    real energy = prefactor*qq * (1.0 / sqrt(distSqr) - B1_half*distSqr -crf);     
                    return energy;
                }

                bool _computeForce(Real3D& force, const Particle &p1,
                         const Particle &p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real r2 = dist.sqr();
                    if (r2>rc2) return true;
                    real r = sqrt(r2);
                    real qq = p1.q()*p2.q();
                    if (checkTIpair(p1.id(),p2.id())) {
                      qq *= complLambdaTI; //(1-lambda)*qAi*qAj
                    }
                    real ffactor = prefactor*qq* (1.0/(r*r2) + B1);
                    force = dist * ffactor;
                    return true;

                }

                real _computeEnergyDeriv(const Particle& p1, const Particle& p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real distSqr = dist.sqr();
                    if (distSqr>rc2) return 0.0;
                    if (checkTIpair(p1.id(),p2.id())) {
                        real qq = p1.q()*p2.q();
                        real dhdl = -1.0 * prefactor * qq * (1.0 / sqrt(distSqr) - B1_half*distSqr -crf);
                        return dhdl;
                    } else {
                        return 0.0;
                    }
                }
                
                real _computeEnergySqrRaw(real distSqr) const {
                        cout << "_computeEnergySqrRaw not possible for reaction field, no particle information" << endl;
                        exit(0);
                        return 0;
                }
                bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
                        cout << "_computeEnergySqrRaw not possible for reaction field, no particle information" << endl;
                        exit(0);
                        return false;
                }
        };


        // provide pickle support
        struct ReactionFieldGeneralizedTI_pickle : boost::python::pickle_suite
        {
          static
          boost::python::tuple
          getinitargs(ReactionFieldGeneralizedTI const& pot)
          {
              real _prefactor = pot.getPrefactor();
              real _kappa     = pot.getKappa();
              real _eps1      = pot.getEpsilon1();
              real _eps2      = pot.getEpsilon2();
              real rc         = pot.getCutoff();
              real _lambda    = pot.getLambdaTI();
              real _annihil   = pot.getAnnihilate();
              return boost::python::make_tuple(_prefactor, _kappa, _eps1, _eps2, rc, _lambda, _annihil);
          }
        };

    }
}

#endif
