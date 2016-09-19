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
#ifndef _INTERACTION_LENNARDJONESSOFTCORETI_HPP
#define _INTERACTION_LENNARDJONESSOFTCORETI_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "Potential.hpp"
#include <set>

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the softcore Lennard Jones potential
        in a Thermodynamic Integration simulation
    **/
    class LennardJonesSoftcoreTI : public PotentialTemplate< LennardJonesSoftcoreTI > {
    private:
      real epsilonA, epsilonB;
      real ff1A, ff2A;
      real ff1B, ff2B;
      bool annihilate; //if true, atoms in pidsTI are annihilated, otherwise they are decoupled
      std::set<longint> pidsTI; //PIDs of particles whose LJ interaction is zero in TI state B
      real lambdaTI;
      real complLambdaTI; //1-lambdaTI
      real alphaSC;
      real sigmaSC_A, sigmaSC_B;
      real powerSC;
      real sigmaSC_A6; //pow(sigmaSC_A,6);
      real sigmaSC_B6; //pow(sigmaSC_B,6);
      real lambdaTI_powerSC; //pow(lambdaTI,powerSC);
      real lambdaTI_powerSCm1; //pow(lambdaTI,powerSC-1.0);
      real compllambdaTI_powerSC; //pow(1.0-lambdaTI,powerSC);
      real compllambdaTI_powerSCm1; //pow(1.0-lambdaTI,powerSC-1.0);
      real alpha_sigmaA6_lambdaP; //alphaSC * sigmaSC_A6 * lambdaTI_powerSC;
      real alpha_sigmaB6_compllambdaP; //alphaSC * sigmaSC_B6 * compllambdaTI_powerSC;
      real powerSC_alphaSC_inv6; //powerSC*alphaSC/6.0;

      void preset() {
        complLambdaTI = 1.0 - lambdaTI;

        sigmaSC_A6 = pow(sigmaSC_A,6);
        sigmaSC_B6 = pow(sigmaSC_B,6);
        lambdaTI_powerSC = pow(lambdaTI,powerSC);
        lambdaTI_powerSCm1 = pow(lambdaTI,powerSC-1.0);
        compllambdaTI_powerSC = pow(1.0-lambdaTI,powerSC);
        compllambdaTI_powerSCm1 = pow(1.0-lambdaTI,powerSC-1.0);
        alpha_sigmaA6_lambdaP = alphaSC * sigmaSC_A6 * lambdaTI_powerSC;
        alpha_sigmaB6_compllambdaP = alphaSC * sigmaSC_B6 * compllambdaTI_powerSC;
        powerSC_alphaSC_inv6 = powerSC*alphaSC/6.0;

        real sig2 = sigmaSC_A * sigmaSC_A; 
        real sig6 = sig2 * sig2 * sig2; 
        ff1A = 48.0 * epsilonA * sig6 * sig6; 
        ff2A = 24.0 * epsilonA * sig6; 

        sig2 = sigmaSC_B * sigmaSC_B; 
        sig6 = sig2 * sig2 * sig2; 
        ff1B = 48.0 * epsilonB * sig6 * sig6; 
        ff2B = 24.0 * epsilonB * sig6; 
      }

    public:
      static void registerPython();

      LennardJonesSoftcoreTI()
	: epsilonA(0.0), sigmaSC_A(0.0), epsilonB(0.0), sigmaSC_B(0.0), alphaSC(0.0), powerSC(0.0),
          lambdaTI(0.0), annihilate(0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      LennardJonesSoftcoreTI(real _epsilonA, real _sigmaA, 
                   real _epsilonB, real _sigmaB,
                   real _alpha, real _p,
		   real _cutoff, real _lambda, bool _annihilate) 
	: epsilonA(_epsilonA), sigmaSC_A(_sigmaA), 
          epsilonB(_epsilonB), sigmaSC_B(_sigmaB),
          alphaSC(_alpha), powerSC(_p),
          lambdaTI(_lambda), annihilate(_annihilate) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
      }

      virtual ~LennardJonesSoftcoreTI() {};

      real getLambdaTI() const { return lambdaTI; }
      real getAlphaSC() const { return alphaSC; }
      real getPowerSC() const { return powerSC; }

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
        if (distSqr>cutoffSqr) return 0.0;

        real frac2, frac6;
        real energyA, energyB;
        real rA, rB, rA5, rB5, rA2, rB2;
        real r6, r5;

        if (checkTIpair(p1.id(),p2.id())) { //use softcore and lambdaTI
          r6 = distSqr*distSqr*distSqr;

          rA = pow(alpha_sigmaA6_lambdaP + r6,1.0/6.0);
          rA2 = rA*rA; 
          frac2 = sigmaSC_A*sigmaSC_A / rA2;
          frac6 = frac2 * frac2 * frac2;
          energyA = 4.0 * epsilonA * (frac6 * frac6 - frac6);

          rB = pow(alpha_sigmaB6_compllambdaP + r6,1.0/6.0);
          rB2 = rB*rB;
          frac2 = sigmaSC_B*sigmaSC_B / rB2;
          frac6 = frac2 * frac2 * frac2;
          energyB = 4.0 * epsilonB * (frac6 * frac6 - frac6);

          real energy = complLambdaTI*energyA + lambdaTI*energyB;
          return energy;
        } else {
          //hardcore, epsilon and sigma used are from state A
          frac2 = sigmaSC_A*sigmaSC_A / distSqr; 
          frac6 = frac2 * frac2 * frac2;
          real energy = 4.0 * epsilonA * (frac6 * frac6 - frac6);
          return energy;
        }
      }

      bool _computeForce(Real3D& force, const Particle &p1,
                         const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        real distSqr = dist.sqr();
        if (distSqr>cutoffSqr) return true;

        real frac2, frac6;

        if (checkTIpair(p1.id(),p2.id())) { 

          real forceA, forceB; //this is dV/dr not dV/dr * 1/r
          real rA, rB, rA5, rB5, rA2, rB2;
          real r6, r5;

          r6 = distSqr*distSqr*distSqr;
          r5 = distSqr*distSqr*sqrt(distSqr);

          rA = pow(alpha_sigmaA6_lambdaP + r6,1.0/6.0);
          rA2 = rA*rA; 
          rA5 = rA2*rA2*rA;
          frac2 = 1.0 / rA2;
          frac6 = frac2 * frac2 * frac2;
          forceA = frac6 * (ff1A * frac6 - ff2A) * frac2 * rA; 

          rB = pow(alpha_sigmaB6_compllambdaP + r6,1.0/6.0);
          rB2 = rB*rB; 
          rB5 = rB2*rB2*rB;
          frac2 = 1.0 / rB2;
          frac6 = frac2 * frac2 * frac2;
          forceB = frac6 * (ff1B * frac6 - ff2B) * frac2 * rB;

          real ffactor = complLambdaTI*r5/rA5*forceA + lambdaTI*r5/rB5*forceB; 
          //ffactor = dV/dr * 1/r; force = dV/dr
          ffactor /= sqrt(distSqr);
          force = dist * ffactor ;
          return true;
        } else {
          real frac2 = 1.0 / distSqr;
          real frac6 = frac2 * frac2 * frac2;
          real ffactor = frac6 * (ff1A * frac6 - ff2A) * frac2; //dV/dr * 1/r
          force = dist * ffactor;
          return true;
        }
      }

      real _computeEnergyDeriv(const Particle& p1, const Particle& p2) const {          
        Real3D dist = p1.position() - p2.position();
        real distSqr = dist.sqr();
        if (distSqr>cutoffSqr) return 0.0;
        if (checkTIpair(p1.id(),p2.id())) {

          real frac2, frac6;
          real energyA, energyB;

          real rA, rB, rA5, rB5, rA2, rB2;
          real r6, r5;

          r6 = distSqr*distSqr*distSqr;

          rA = pow(alpha_sigmaA6_lambdaP + r6,1.0/6.0);
          rA2 = rA*rA; 
          frac2 = sigmaSC_A*sigmaSC_A / rA2;
          frac6 = frac2 * frac2 * frac2;
          energyA = 4.0 * epsilonA * (frac6 * frac6 - frac6);

          rB = pow(alpha_sigmaB6_compllambdaP + r6,1.0/6.0);
          rB2 = rB*rB;
          frac2 = sigmaSC_B*sigmaSC_B / rB2;
          frac6 = frac2 * frac2 * frac2;
          energyB = 4.0 * epsilonB * (frac6 * frac6 - frac6);

          real forceA, forceB; //dV/dr not dV/dr * 1/r

          r5 = distSqr*distSqr*sqrt(distSqr);

          rA5 = rA2*rA2*rA;
          frac2 = 1.0 / rA2;
          frac6 = frac2 * frac2 * frac2;
          forceA = frac6 * (ff1A * frac6 - ff2A) * frac2 * rA; 

          rB5 = rB2*rB2*rB;
          frac2 = 1.0 / rB2;
          frac6 = frac2 * frac2 * frac2;
          forceB = frac6 * (ff1B * frac6 - ff2B) * frac2 * rB;

          real energyDeriv = energyB - energyA + powerSC_alphaSC_inv6 * (
           lambdaTI * forceB / rB5 * sigmaSC_B6 * compllambdaTI_powerSCm1
           - complLambdaTI * forceA / rA5 * sigmaSC_A6 * lambdaTI_powerSCm1);

          return energyDeriv;
        } else {
          return 0.0;
        }
      }

      real _computeEnergySqrRaw(real distSqr) const {
              std::cout << "_computeEnergySqrRaw not implemented" << std::endl;
              exit(0);
              return 0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
              std::cout << "_computeEnergySqrRaw not implemented" << std::endl;
              exit(0);
              return false;
      }
    };

  }
}

#endif
