// ESPP_CLASS
#ifndef _INTERACTION_REACTIONFIELDGENERALIZED_HPP
#define _INTERACTION_REACTIONFIELDGENERALIZED_HPP

#include "Potential.hpp"
//#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
    namespace interaction {
        /*
         * This class provides methods to compute forces and energies of
         * the generalized reaction field
         */
        class ReactionFieldGeneralized: public
        PotentialTemplate<ReactionFieldGeneralized> {
            private:
                real qq; // p1.q() * p2.q()
                real kappa;
                real epsilon1, epsilon2;
                real rc, rc3; // cutoff, cutoff^3
                real B1, B0;
                real prefactor;

                void initialize() {
                    real krc = kappa*rc;
                    real tmp1 = (epsilon1 - 4*epsilon2) *
                            (1 + krc) - 2*epsilon2 * krc*krc;
                    real tmp2 = (epsilon1 + 2*epsilon2) *
                            (1 + krc) + epsilon2 * krc*krc;
                    rc3 = pow(rc,3);
                    B1 = tmp1/tmp2;
                    B1 = (1+B1) / rc3;
                    B0 = epsilon1 - 2*epsilon2 * (1+krc);
                    B0 = (1+B0) / rc;
                    prefactor = 176.788 * 3.9176;
                    qq = qq*prefactor / epsilon1;
                }

            public:
                static void registerPython();

                ReactionFieldGeneralized()
                :qq(0.0), kappa(0.0),
                 epsilon1(1.0), epsilon2(80.0),
                 rc(1.0){
                    setShift(0.0);
                    setCutoff(infinity);
                    initialize();
                }

                ReactionFieldGeneralized(
                        real _qq,
                        real _kappa,
                        real _eps1,
                        real _eps2,
                        real _cutoff,
                        real _shift)
                : qq(_qq), kappa(_kappa),
                  epsilon1(_eps1), epsilon2(_eps2),
                  rc(_cutoff) {
                    setShift(_shift);
                    setCutoff(_cutoff);

                    initialize();
                }

                ReactionFieldGeneralized(
                        real _qq,
                        real _kappa,
                        real _eps1,
                        real _eps2,
                        real _cutoff)
                : qq(_qq), kappa(_kappa),
                  epsilon1(_eps1), epsilon2(_eps2),
                  rc(_cutoff) {
                    autoShift = false;
                    setCutoff(_cutoff);
                    setAutoShift();

                    initialize();
                }

                void setCutoff(real _cutoff) {
                    cutoff = _cutoff;
                    cutoffSqr = cutoff*cutoff;
                    updateAutoShift();
                }

                // Setter and getter
                void setQQ(real _qq) {
                    qq = _qq;
                    updateAutoShift();
                }
                real getQQ() const { return qq; }


                /*
                real _computeEnergy(const Particle& p1, const Particle& p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real distSqr = dist.sqr();
                    real qq = p1.q()*p2.q() / epsilon1;
                    real energy = qq * _computeEnergySqrRaw(distSqr);

                    return energy;
                }*/


                real _computeEnergySqrRaw(real distSqr) const {
                    return qq * (1.0 / sqrt(distSqr) + B0);
                }


                /*
                bool _computeForce(Real3D& force, const Particle &p1,
                         const Particle &p2) const {
                    Real3D dist = p1.position() - p2.position();
                    real r2 = dist.sqr();
                    real r = sqrt(r2);

                    real qq = p1.q()*p2.q() / epsilon1;
                    real ffactor = qq * (1/(r*r2) + (1+B1)/rc3);
                    force = dist * ffactor;

                    return true;
                }*/

                bool _computeForceRaw(Real3D& force,
                        const Real3D& dist, real distSqr) const {

                    real r = sqrt(distSqr);
                    real ffactor = qq * ((1/(r*distSqr)) + B1);
                    force = dist * ffactor;

                    return true;
                }
        };
    }
}

#endif
