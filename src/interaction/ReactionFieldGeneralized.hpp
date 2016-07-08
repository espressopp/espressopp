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
#ifndef _INTERACTION_REACTIONFIELDGENERALIZED_HPP
#define _INTERACTION_REACTIONFIELDGENERALIZED_HPP

#include "Potential.hpp"
//#include "FixedPairListInteractionTemplate.hpp"

using namespace std;

namespace espressopp {
    namespace interaction {
        /*
         * This class provides methods to compute forces and energies of
         * the generalized reaction field
         */
        class ReactionFieldGeneralized: public
        PotentialTemplate<ReactionFieldGeneralized> {
            private:
                //real qq; // p1.q() * p2.q()
                real kappa;
                real epsilon1, epsilon2;
                real rc, rc3; // cutoff, cutoff^3
                real rc2;
                real B1, B0, B1_half;
                real prefactor;
                real crf; /* const to make potential zero at cutoff, corresponding
                 to the shift, which we cannot use because it depends on qq*/

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
                    //B0 = (epsilon1 - 2.0*epsilon2 * (1.0+krc))/(epsilon2*(1.0+krc));
                    crf = 3*epsilon2/(rc*(2*epsilon2+1.0));
                    //B0 = (1.0+B0) / rc;
                    prefactor/=epsilon1;
                }

            public:
                static void registerPython();

                ReactionFieldGeneralized()
                : prefactor(0.0), kappa(0.0),
                 epsilon1(1.0), epsilon2(80.0),
                 rc(1.0){
                    setShift(0.0);
                    setCutoff(infinity);
                    initialize();
                }

                ReactionFieldGeneralized(   
                        real _prefactor,
                        real _kappa,
                        real _eps1,
                        real _eps2,
                        real _cutoff,
                        real _shift)
                : prefactor(_prefactor),kappa(_kappa),
                  epsilon1(_eps1), epsilon2(_eps2),
                  rc(_cutoff) {
                    setShift(_shift);//setShift(_shift);
                    setCutoff(_cutoff);

                    initialize();
                }

                ReactionFieldGeneralized(
                        real _prefactor,
                        real _kappa,
                        real _eps1,
                        real _eps2,
                        real _cutoff)
                : prefactor(_prefactor), kappa(_kappa),
                  epsilon1(_eps1), epsilon2(_eps2),
                  rc(_cutoff) {
                    autoShift = false;
                    setCutoff(_cutoff);
                    /*Note: AutoShift cannot be used here since the shift
                     *  has to depend on the product of charges */
                    //setAutoShift();                 
                    initialize();
                }

                real getKappa() const { return kappa; }
                real getEpsilon1() const { return epsilon1; }
                real getEpsilon2() const { return epsilon2; }

                void setCutoff(real _cutoff) {
                    cutoff = _cutoff;
                    cutoffSqr = cutoff*cutoff;
                    //updateAutoShift();
                }

                // Setter and getter
                void setPrefactor(real _prefactor) {
                    prefactor = _prefactor;
                    //updateAutoShift();
                }
                real getPrefactor() const { return prefactor; }


                
                real _computeEnergy(const Particle& p1, const Particle& p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real distSqr = dist.sqr();
                    //std::cout << "qq" << distSqr << " " << p1.id() << " " << p2.id() << "\n";
                    if (distSqr>rc2) return 0.0;
                    real qq = p1.q()*p2.q();
                    /* Note: this implementation counts minus integral of the force as energy
                         which is not the same as the full electrostatic energy
                         See the original paper Tironi et al, J. Chem. Phys. 102 , 5451 (1995) 
                         for details*/

                    
                    
                    real energy = prefactor*qq * (1.0 / sqrt(distSqr) - B1_half*distSqr -crf);     
                    return energy;

                    // FORCE CAPPING HACK (was temporarily used for some ideal gas test simulations)
                    /*real caprad = 0.1;
                    real capradSqr = caprad * caprad;
                    
                    if (distSqr > capradSqr) {
                        real energy = prefactor*qq * (1.0 / sqrt(distSqr) - B1_half*distSqr -crf);     
                        return energy;                
                    }
                    else{
                        real energy = prefactor*qq * (1.0 / sqrt(capradSqr) - B1_half*capradSqr -crf);
                        real forcepart = prefactor*qq* (1.0/(caprad*capradSqr) + B1) *caprad;
                        real out = energy + forcepart*(caprad-sqrt(distSqr));
                        return out;                                        
                    }*/
             
                }


                /*real _computeEnergySqrRaw(real distSqr) const {
                    return qq * (1.0 / sqrt(distSqr) + B0);
                }*/


                
                bool _computeForce(Real3D& force, const Particle &p1,
                         const Particle &p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real r2 = dist.sqr();
                    if (r2>rc2) return true;
                    
                    
                                      
                    real r = sqrt(r2);
                    real qq = p1.q()*p2.q();
                    real ffactor = prefactor*qq* (1.0/(r*r2) + B1);
                    force = dist * ffactor;
                    return true;
                               
                    // FORCE CAPPING HACK (was temporarily used for some ideal gas test simulations)
                    /*real caprad = 0.1;
                    real capradSqr = caprad * caprad;
                    real r = sqrt(r2);
                    real qq = p1.q()*p2.q();
                    
                    if (r2 > capradSqr) {
                        real ffactor = prefactor*qq* (1.0/(r*r2) + B1);
                        force = dist * ffactor;
                        return true;                
                    }
                    else{
                        real ffactor = prefactor*qq* (1.0/(caprad*capradSqr) + B1);
                        force = dist * ffactor * (caprad/r);
                        //std::cout << "ReactionField, capped Force: " << sqrt(dist.sqr()) * ffactor * (caprad/r) << "\n"; 0.1 LEADS TO 2492.93
                        return true;                                         
                    }*/

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
                /*bool _computeForceRaw(Real3D& force,
                        const Real3D& dist, real distSqr) const {

                    real r = sqrt(distSqr);
                    real ffactor = qq * ((1/(r*distSqr)) + B1);
                    force = dist * ffactor;

                    return true;
                }*/
        };

        // provide pickle support
        struct ReactionFieldGeneralized_pickle : boost::python::pickle_suite
        {
          static
          boost::python::tuple
          getinitargs(ReactionFieldGeneralized const& pot)
          {
              real _prefactor = pot.getPrefactor();
              real _kappa     = pot.getKappa();
              real _eps1      = pot.getEpsilon1();
              real _eps2      = pot.getEpsilon2();
              real rc         = pot.getCutoff();
              real sh         = pot.getShift();
              return boost::python::make_tuple(_prefactor, _kappa, _eps1, _eps2, rc, sh);
          }
        };

    }
}

#endif
