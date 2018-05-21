/*
  Copyright (C) 2017,2018
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

#include "python.hpp"
#include "PIAdressIntegrator.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "VerletListAdress.hpp"

namespace espressopp {
  using namespace std;
  namespace integrator {
    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    PIAdressIntegrator::PIAdressIntegrator(shared_ptr< System > system, shared_ptr<VerletListAdress> _verletList)
      : MDIntegrator(system), verletList(_verletList)
    {
      // Initialize variables
      resortFlag = true;
      maxDist    = 0.0;
      verletlistBuilds = 0;

      sStep = 1;
      mStep = 1;
      dt2 = dt*sStep;
      dt3 = dt*mStep*sStep;
      temperature = 0.0;
      gamma = 0.0;
      ntrotter = 1;
      CMDparameter = 1.0;
      PILElambda = 0.5;
      speedup = false;
      constkinmass = false;
      realkinmass = true;
      KTI = false;
      centroidthermostat = true;
      PILE = false;
      clmassmultiplier = 100.0;

      setOmegaSquared();

      // AdResS geometry
      dhy = verletList->getHy();
      pidhy2 = M_PI/(dhy * 2.0);
      dex = verletList->getEx();
      dex2 = dex * dex;
      dexdhy = dex + verletList->getHy();
      dexdhy2 = dexdhy * dexdhy;

      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }
      rng = system->rng;
    }


    PIAdressIntegrator::~PIAdressIntegrator() {}


    void PIAdressIntegrator::run(int nsteps) {
      if(Eigenvalues[0] > 0.00000000001) {
        throw std::runtime_error("Eigenvalues don't start with zero!");
      }

      if(PILE == true && realkinmass == false) {
        throw std::invalid_argument("Using the path integral Langevin scheme (PILE) only makes sense when using real masses for kinetic masses!");
      }

      System& system = getSystemRef();
      storage::Storage& storage = *system.storage;
      real skinHalf = 0.5 * system.getSkin();

      // Set the centroids, weights, and masses
      initializeSetup();

      if (resortFlag) {
        storage.decompose();
        verletlistBuilds += 1;
        maxDist = 0.0;
        resortFlag = false;
      }

      // Update mode positions
      transPos2();

      // Update mode momenta
      transMom2();

      // outer loop (interatomic nonbonded forces)
      for (int i = 0; i < nsteps; i++) {
        if(i==0) {
          updateForces(3);
          integrateV1(3, false);
        }

        // middle loop (interatomic bonded forces)
        for (int j = 0; j < mStep; j++) {
          if(j==0) {
            updateForces(2);
            integrateV1(2, false);
          }

          // inner loop (intraatomic mode forces, thermostat and position and momenta propagation)
          for (int k= 0; k < sStep; k++) {
            if(k==0) {
              updateForces(1);
            }
            integrateV2();

            integrateModePos();
            OUintegrate();
            integrateModePos();

            updateForces(1);
            integrateV2();
          }

          // Update real positions (inner loop works only in mode space)
          transPos1();

          updateForces(2);
          if(j==(mStep-1)) {
            integrateV1(2,false);
          }
          else {
            integrateV1(2,true);
          }
        }

        // If necessary, rebuild the Verlet list
        if (maxDist > skinHalf) resortFlag = true;
        if (resortFlag) {
          storage.decompose();
          verletlistBuilds += 1;
          // Update mode positions after rebuild
          transPos2();
          maxDist  = 0.0;
          resortFlag = false;
        }

        updateForces(3);
        if(i==(nsteps-1)) {
          integrateV1(3,false);
        }
        else {
          integrateV1(3,true);
        }
        step++;
      }

      // Update real velocities at the end
      transMom1();

    }


    using namespace boost::python;


    void PIAdressIntegrator::initializeSetup() {
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);

        if (it3 != fixedtupleList->end()) {

          std::vector<Particle*> atList;
          atList = it3->second;

          Real3D cmp(0.0, 0.0, 0.0);
          Real3D cmv(0.0, 0.0, 0.0);
          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            cmp += at.position();
            cmv += at.velocity();
          }
          cmv /= ntrotter;
          cmp /= ntrotter;

          // update (overwrite) the position and velocity of the VP
          vp.position() = cmp;
          vp.velocity() = cmv;

          if (KTI == false) {
            std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
            Real3D pa = **it2; // position of adress particle
            Real3D d1(0.0, 0.0, 0.0);
            verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
            real min1sq;
            if (verletList->getAdrRegionType()) { // spherical adress region
              min1sq = d1.sqr();
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                pa = **it2;
                verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                real distsq1 = d1.sqr();
                if (distsq1 < min1sq) min1sq = distsq1;
              }
            }
            else { //slab-type adress region
              min1sq = d1[0]*d1[0];
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                pa = **it2; // position of adress particle
                verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                real distsq1 = d1[0]*d1[0];
                if (distsq1 < min1sq) min1sq = distsq1;
              }
            }
            real w = weight(min1sq);
            vp.lambda() = w;
            real wDeriv = weightderivative(sqrt(min1sq));
            vp.lambdaDeriv() = wDeriv;
            vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
          }

        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
    }


    void PIAdressIntegrator::integrateV1(int t, bool doubletime) {
      real half_dt = 0.0;
      if (t==2) {
        if(doubletime) {
          half_dt = dt2;
        }
        else {
          half_dt = 0.5 * dt2;
        }
      }
      else if (t==3) {
        if(doubletime) {
          half_dt = dt3;
        }
        else {
          half_dt = 0.5 * dt3;
        }
      }
      else {
        throw std::runtime_error("integrateV1 routine in PIAdressIntegrator received wrong integer.");
      }

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;

          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            if((at.pib() != 1) && (speedup == true) && (vp.lambda() < 0.000000001)) {
              continue;
            } else {
              at.modemom() += half_dt * at.forcem();
              if(at.pib() == 1)
              {
                vp.velocity() = (1.0/sqrt(ntrotter)) * at.modemom()/(vp.mass());
              }
            }
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }

    }


    void PIAdressIntegrator::integrateV2() {
      real half_dt = 0.5 * dt;
      real half_dt4 = 0.25 * dt;

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;

          if (constkinmass == false) { // kinetic mass also changes - hence, second decomposition level required

            // First half_dt4 loop without centroid mode
            if((speedup == false) || (vp.lambda() > 0.000000001)) { // do not propagate the higher modes in the CL region when using speedup
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                   it2 != atList.end(); ++it2) {
                Particle &at = **it2;
                if(at.pib() == 1) { // exclude centroid
                  continue;
                }
                else if(at.pib() > 1 && at.pib() <= ntrotter) {
                  at.modemom() += half_dt4 * at.forcem();
                }
                else {
                  throw std::runtime_error("at.pib() outside of trotter range in integrateV2 function (PIAdressIntegrator)!");
                }
              }
            }

            // half_dt2, only centroid mode
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;
              if(at.pib() == 1) { // only centroid

                // If necessary, calculate drift force
                if(vp.lambda()<1.0 && vp.lambda()>0.0) {

                  // Get the inner factor
                  real xi = 0.0;
                  for (std::vector<Particle*>::iterator it5 = atList.begin();
                       it5 != atList.end(); ++it5) {
                    Particle &at2 = **it5;
                    if(at2.pib() != 1) {
                      if(realkinmass == false) {
                        xi += at2.modemom().sqr()/(Eigenvalues[at2.pib()-1]);
                      }
                      else {
                        xi += at2.modemom().sqr();
                      }
                    }
                  }

                  std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                  Real3D pa = **it2;
                  Real3D mindriftforce(0.0, 0.0, 0.0);
                  verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                  real min1sq = 0.0;
                  if(verletList->getAdrRegionType()) {
                    min1sq = mindriftforce.sqr();
                  }
                  else {
                    min1sq = mindriftforce[0]*mindriftforce[0];
                  }
                  ++it2;
                  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                    pa = **it2;
                    Real3D driftforce(0.0, 0.0, 0.0);
                    verletList->getSystem()->bc->getMinimumImageVector(driftforce, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                    real distsq1 = 0.0;
                    if(verletList->getAdrRegionType()) {
                      distsq1 = driftforce.sqr();
                    }
                    else {
                      distsq1 = driftforce[0]*driftforce[0];
                    }
                    if (distsq1 < min1sq) {
                      min1sq = distsq1;
                      mindriftforce = driftforce;
                    }
                  }
                  min1sq = sqrt(min1sq);
                  if(verletList->getAdrRegionType()) {
                    mindriftforce = (1.0/min1sq)*mindriftforce;
                    mindriftforce *= 0.5*(clmassmultiplier-1.0)*vp.mass()*sqrt(ntrotter)*xi*vp.lambdaDeriv()/(vp.varmass()*vp.varmass()*CMDparameter);
                    at.modemom() += half_dt * (at.forcem() - mindriftforce);
                  }
                  else {
                    real mindriftforceX = (1.0/min1sq)*mindriftforce[0];
                    mindriftforceX *= 0.5*(clmassmultiplier-1.0)*vp.mass()*sqrt(ntrotter)*xi*vp.lambdaDeriv()/(vp.varmass()*vp.varmass()*CMDparameter);
                    Real3D driftforceadd(mindriftforceX,0.0,0.0);
                    at.modemom() += half_dt * (at.forcem() - driftforceadd);
                  }
                }
                else {
                  at.modemom() += half_dt * at.forcem();
                }
                vp.velocity() = sqrt(ntrotter) * at.modemom()/(vp.mass());
                break; // Once we propagated the centroid, we can leave the loop

              }
              else if(at.pib() > 1 && at.pib() <= ntrotter) {
                continue;
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in integrateV2 function (PIAdressIntegrator)!");
              }
            }

            // Second half_dt4 loop without centroid mode
            if((speedup == false) || (vp.lambda() > 0.000000001)) { // do not propagate the higher modes in the CL region when using speedup

              for (std::vector<Particle*>::iterator it2 = atList.begin();
                   it2 != atList.end(); ++it2) {
                Particle &at = **it2;

                if(at.pib() == 1) { // exclude centroid
                  continue;
                }
                else if(at.pib() > 1 && at.pib() <= ntrotter) {
                  at.modemom() += half_dt4 * at.forcem();
                }
                else {
                  throw std::runtime_error("at.pib() outside of trotter range in integrateV2 function (PIAdressIntegrator)!");
                }
              }
            }
          }
          else { // constant kinetic mass - hence, no second decomposition level required
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;
              if((at.pib() != 1) && (speedup == true) && (vp.lambda() < 0.000000001)) {
                continue;
              } else {
                at.modemom() += half_dt * at.forcem();
                if(at.pib() == 1)
                {
                  vp.velocity() = sqrt(ntrotter) * at.modemom()/(vp.mass());
                }
              }
            }
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }

    }


    void PIAdressIntegrator::integrateModePos() {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells(); //!!!!
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;

        if (constkinmass == false) { // kinetic mass also changes - hence, second decomposition level required

          real half_dt4 = 0.25 * ntrotter *  dt / (vp.mass());
          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);
          if (it3 != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it3->second;

            // First half_dt4 loop with only centroid mode
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;
              if(at.pib() == 1) { // only centroid
                at.modepos() += half_dt4 * at.modemom();

                if (KTI == false) {
                  // Update resolution and variable masses
                  std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                  Real3D pa = **it2;
                  Real3D d1(0.0, 0.0, 0.0);
                  real min1sq;
                  verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                  if (verletList->getAdrRegionType()) { // spherical adress region
                    min1sq = d1.sqr();
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = d1.sqr();
                      if (distsq1 < min1sq) min1sq = distsq1;
                    }
                  }
                  else { //slab-type adress region
                    min1sq = d1[0]*d1[0];
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = d1[0]*d1[0];
                      if (distsq1 < min1sq) min1sq = distsq1;
                    }
                  }
                  real w = weight(min1sq);
                  vp.lambda() = w;
                  real wDeriv = weightderivative(min1sq);
                  vp.lambdaDeriv() = wDeriv;
                  vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
                }

                break; // Once we propagated the centroid, we can leave the loop
              } else {
                continue;
              }
            }

            // half_dt2 loop without centroid mode
            real half_dt = 0.5 * ntrotter * dt / (CMDparameter * vp.varmass());
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;
              if((at.pib() == 1) || ((speedup == true) && (vp.lambda() < 0.000000001) && (at.pib() != 1))) { // do not propagate the higher modes in the CL region when using speedup, also exclude centroid
                continue;
              }
              else {
                if(realkinmass == false) {
                  at.modepos() += half_dt * at.modemom() / (Eigenvalues[at.pib()-1]);
                }
                else {
                  at.modepos() += half_dt * at.modemom();
                }
              }
            }

            // Second half_dt4 loop with only centroid mode
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;

              if(at.pib() == 1) { // only centroid
                at.modepos() += half_dt4 * at.modemom();

                if (KTI == false) {
                  // Update resolution and variable masses
                  std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                  Real3D pa = **it2;
                  Real3D d1(0.0, 0.0, 0.0);
                  real min1sq;
                  verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                  if (verletList->getAdrRegionType()) { // spherical adress region
                    min1sq = d1.sqr();
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = d1.sqr();
                      if (distsq1 < min1sq) min1sq = distsq1;
                    }
                  }
                  else { //slab-type adress region
                    min1sq = d1[0]*d1[0];
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = d1[0]*d1[0];
                      if (distsq1 < min1sq) min1sq = distsq1;
                    }
                  }
                  real w = weight(min1sq);
                  vp.lambda() = w;
                  real wDeriv = weightderivative(min1sq);
                  vp.lambdaDeriv() = wDeriv;
                  vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
                }

                break; // Once we propagated the centroid, we can leave the loop
              } else {
                continue;
              }
            }
          }
          else {
            std::stringstream ss;
            ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
            throw std::runtime_error(ss.str());
          }

        }
        else { // constant kinetic mass - hence, no second decomposition level required

          real half_dt = 0.5 * ntrotter * dt / (vp.mass());
          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);
          if (it3 != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it3->second;

            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;

              if(at.pib() == 1) {
                at.modepos() += half_dt * at.modemom();

                if (KTI == false) {
                  // Update resolution and variable masses
                  std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                  Real3D pa = **it2;
                  Real3D d1(0.0, 0.0, 0.0);
                  real min1sq;
                  verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);

                  if (verletList->getAdrRegionType()) { // spherical adress region
                    min1sq = d1.sqr();
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = d1.sqr();
                      if (distsq1 < min1sq) min1sq = distsq1;
                    }
                  }
                  else { //slab-type adress region
                    min1sq = d1[0]*d1[0];
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = d1[0]*d1[0];
                      if (distsq1 < min1sq) min1sq = distsq1;
                    }
                  }
                  real w = weight(min1sq);
                  vp.lambda() = w;
                  real wDeriv = weightderivative(min1sq);
                  vp.lambdaDeriv() = wDeriv;
                  vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
                }

              } else if(at.pib() > 1 && at.pib() <= ntrotter) {
                if((speedup == true) && (vp.lambda() < 0.000000001)) {
                  continue;
                } else {
                  if(realkinmass == false) {
                    at.modepos() += half_dt * at.modemom() / (CMDparameter * Eigenvalues[at.pib()-1]);
                  }
                  else {
                    at.modepos() += half_dt * at.modemom() / CMDparameter;
                  }
                }
              } else {
                throw std::runtime_error("at.pib() outside of trotter range in integrateModePos function (PIAdressIntegrator)!");
              }
            }
          }
          else {
            std::stringstream ss;
            ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
            throw std::runtime_error(ss.str());
          }

        }
      }
    }


    void PIAdressIntegrator::OUintegrate() { // this is the Ornstein-Uhlenbeck process, this is, the Langevin thermostat
      real prefac1 = exp(-gamma*dt);
      real prefac2 = sqrt(12.0*temperature*( 1.0-exp(-2.0*gamma*dt) ) / ntrotter);
      // careful with 12... it's because the uniform distribution needs to have the same width as the Gaussian one.
      // also note that kb is inside temperature, we use gromacs units

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;

          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            if(at.pib() == 1) {
              if(centroidthermostat == true || vp.lambda() < 1.0) {
                Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
                at.modemom() = prefac1 * at.modemom() + prefac2*sqrt(vp.mass()) * ranval;
                vp.velocity() = sqrt(ntrotter) * at.modemom()/(vp.mass());
              }
            }
            else if(at.pib() > 1 && at.pib() <= ntrotter) {
              if((speedup == false) || (vp.lambda() > 0.000000001)) {
                Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
                if(PILE == false) {
                  if(constkinmass == false) {
                    if(realkinmass == false) {
                      at.modemom() = prefac1 * at.modemom() + (prefac2*sqrt(CMDparameter*vp.varmass()*Eigenvalues[at.pib()-1])) * ranval;
                    }
                    else {
                      at.modemom() = prefac1 * at.modemom() + (prefac2*sqrt(CMDparameter*vp.varmass())) * ranval;
                    }
                  }
                  else {
                    if(realkinmass == false) {
                      at.modemom() = prefac1 * at.modemom() + (prefac2*sqrt(CMDparameter*vp.mass()*Eigenvalues[at.pib()-1])) * ranval;
                    }
                    else {
                      at.modemom() = prefac1 * at.modemom() + (prefac2*sqrt(CMDparameter*vp.mass())) * ranval;
                    }
                  }
                }
                else {
                  real modegamma = 2.0 * PILElambda * sqrt(omega2 * ntrotter) * sqrt(Eigenvalues[at.pib()-1]);
                  real prefac1_PILE = exp(-modegamma*dt);
                  real prefac2_PILE = sqrt(12.0*temperature*( 1.0-exp(-2.0*modegamma*dt) ) / ntrotter );
                  if(constkinmass == false) {
                    at.modemom() = prefac1_PILE * at.modemom() + (prefac2_PILE*sqrt(CMDparameter*vp.varmass())) * ranval;
                  }
                  else {
                    at.modemom() = prefac1_PILE * at.modemom() + (prefac2_PILE*sqrt(CMDparameter*vp.mass())) * ranval;
                  }
                }
              }
            }
            else {
              throw std::runtime_error("at.pib() outside of trotter range in OUintegrate function (PIAdressIntegrator)!");
            }
          }

        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
    }


    void PIAdressIntegrator::transPos1() { // Update the real positions from mode positions
      real maxSqDist = 0.0;
      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if(it2 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it2->second;
          Real3D oldpos = vp.position();

          for (std::vector<Particle*>::iterator it3 = atList.begin();
               it3 != atList.end(); ++it3) {
            Particle &at = **it3;
            Real3D zero(0.0,0.0,0.0);
            at.position()=zero;

            for (std::vector<Particle*>::iterator it5 = atList.begin();
                 it5 != atList.end(); ++it5) {
              Particle &at2 = **it5;
              if(at.pib() <= ntrotter) {
                at.position()+= at2.modepos()*Tvectors[at.pib()-1][at2.pib()-1];
                if((at2.pib() == 1) && (at.pib() == 1)) {
                  vp.position() = (1.0/sqrt(ntrotter)) * at2.modepos();
                  vp.velocity() = sqrt(ntrotter) * at2.modemom()/(vp.mass());
                }
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in transPos1 function (PIAdressIntegrator)!");
              }
            }

          }

          real sqDist = (oldpos-vp.position()).sqr();
          maxSqDist = std::max(maxSqDist, sqDist);
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }

      real maxAllSqDist;
      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
      maxDist+=sqrt(maxAllSqDist);
    }


    void PIAdressIntegrator::transPos2() { // Update the mode positions from real positions
      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if(it2 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it2->second;

          for (std::vector<Particle*>::iterator it3 = atList.begin();
               it3 != atList.end(); ++it3) {
            Particle &at = **it3;
            Real3D zero(0.0,0.0,0.0);
            at.modepos()=zero;

            for (std::vector<Particle*>::iterator it5 = atList.begin();
                 it5 != atList.end(); ++it5) {
              Particle &at2 = **it5;
              if(at.pib() <= ntrotter) {
                at.modepos()+= at2.position()*Eigenvectors[at.pib()-1][at2.pib()-1];
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in transPos2 function (PIAdressIntegrator)!");
              }
            }

          }

        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
    }


    void PIAdressIntegrator::transMom1() { // Update the real velocities from mode momenta
      real maxSqDist = 0.0;

      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if(it2 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it2->second;
          for (std::vector<Particle*>::iterator it3 = atList.begin();
               it3 != atList.end(); ++it3) {
            Particle &at = **it3;
            Real3D zero(0.0,0.0,0.0);
            at.velocity()=zero;

            for (std::vector<Particle*>::iterator it5 = atList.begin();
                 it5 != atList.end(); ++it5) {
              Particle &at2 = **it5;
              if(at.pib() <= ntrotter) {
                if(at2.pib() > 1 && at2.pib() <= ntrotter) {
                  if(constkinmass == false) {
                    if(realkinmass == false) {
                      at.velocity()+= at2.modemom()*Tvectors[at.pib()-1][at2.pib()-1]*ntrotter/((Eigenvalues[at2.pib()-1])*CMDparameter*vp.varmass());
                    }
                    else {
                      at.velocity()+= at2.modemom()*Tvectors[at.pib()-1][at2.pib()-1]*ntrotter/(CMDparameter*vp.varmass());
                    }
                  }
                  else {
                    if(realkinmass == false) {
                      at.velocity()+= at2.modemom()*Tvectors[at.pib()-1][at2.pib()-1]*ntrotter/((Eigenvalues[at2.pib()-1])*CMDparameter*vp.mass());
                    }
                    else {
                      at.velocity()+= at2.modemom()*Tvectors[at.pib()-1][at2.pib()-1]*ntrotter/(CMDparameter*vp.mass());
                    }
                  }
                }
                else if(at2.pib() == 1) {
                  at.velocity()+= at2.modemom()*Tvectors[at.pib()-1][at2.pib()-1]*ntrotter/(vp.mass());
                }
                else {
                  throw std::runtime_error("at.pib() outside of trotter range in transMom1 function (PIAdressIntegrator)!");
                }
                if((at2.pib() == 1) && (at.pib() == 1)) {
                  vp.velocity() = sqrt(ntrotter) * at2.modemom()/(vp.mass());
                }

              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in transMom1 function (PIAdressIntegrator)!");
              }
            }

          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }

      real maxAllSqDist;
      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
      maxDist+=sqrt(maxAllSqDist);
    }


    void PIAdressIntegrator::transMom2() { // Update the mode momenta from real velocities
      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if(it2 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it2->second;
          for (std::vector<Particle*>::iterator it3 = atList.begin();
               it3 != atList.end(); ++it3) {
            Particle &at = **it3;
            Real3D zero(0.0,0.0,0.0);
            at.modemom()=zero;

            for (std::vector<Particle*>::iterator it5 = atList.begin();
                 it5 != atList.end(); ++it5) {
              Particle &at2 = **it5;
              if(at.pib() <= ntrotter) {
                at.modemom()+= at2.velocity()*Eigenvectors[at.pib()-1][at2.pib()-1];
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in transMom2 function (PIAdressIntegrator)!");
              }
            }

            if(at.pib() > 1 && at.pib() <= ntrotter) {
              if(constkinmass == false) {
                if(realkinmass == false) {
                  at.modemom() *= (vp.varmass()*CMDparameter* (Eigenvalues[at.pib()-1])/ntrotter);
                }
                else {
                  at.modemom() *= (vp.varmass()*CMDparameter/ntrotter);
                }
              }
              else {
                if(realkinmass == false) {
                  at.modemom() *= (vp.mass()*CMDparameter* (Eigenvalues[at.pib()-1])/ntrotter);
                }
                else {
                  at.modemom() *= (vp.mass()*CMDparameter/ntrotter);
                }
              }
            }
            else if(at.pib() == 1) {
              at.modemom() *= vp.mass()/ntrotter;
            }
            else {
              throw std::runtime_error("at.pib() outside of trotter range in transMom2 function (PIAdressIntegrator)!");
            }
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
    }


    void PIAdressIntegrator::transForces() { // Update the mode forces from real forces
      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if(it2 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it2->second;
          for (std::vector<Particle*>::iterator it3 = atList.begin();
               it3 != atList.end(); ++it3) {
            Particle &at = **it3;

            for (std::vector<Particle*>::iterator it5 = atList.begin();
                 it5 != atList.end(); ++it5) {
              Particle &at2 = **it5;
              if(at.pib() == 1) {
                at.forcem()+= (1.0/sqrt(ntrotter))* at2.force();
              }
              else if(at.pib() > 1 && at.pib() <= ntrotter) {
                at.forcem()+= at2.force()*Eigenvectors[at.pib()-1][at2.pib()-1];
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in transForces function (PIAdressIntegrator)!");
              }
            }

          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
    }


    void PIAdressIntegrator::calcForcesS() { // Calculate slow forces, this is, all interatomic non-bonded forces
      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;
      for (size_t i = 0; i < srIL.size(); i++) {
        if(srIL[i]->bondType() == Nonbonded) {
          srIL[i]->addForces();
        }
      }
    }


    void PIAdressIntegrator::calcForcesM() { // Calculate medium forces, this is, all interatomic bonded forces (pair and angular)
      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;
      for (size_t i = 0; i < srIL.size(); i++) {
        if(srIL[i]->bondType() == Pair || srIL[i]->bondType() == Angular) {
          srIL[i]->addForces();
        }
      }

      // signal (used for FEC)
      aftCalcF();
    }


    void PIAdressIntegrator::calcForcesF() { // Calculate fast forces, this is, the intra-ring forces between the path integral beads (purely in mode space)
      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if(speedup==true) {

          if(vp.lambda() < 0.000000001) { // In speedup mode we do not integrate the higher modes in the CL region
            continue;
          }
          else {
            if(it2 != fixedtupleList->end()) {
              std::vector<Particle*> atList;
              atList = it2->second;
              for (std::vector<Particle*>::iterator it3 = atList.begin();
                   it3 != atList.end(); ++it3) {
                Particle &at = **it3;
                if(at.pib() > 1 && at.pib() <= ntrotter) {
                  at.forcem() -= at.modepos()*omega2*vp.varmass()*Eigenvalues[at.pib()-1];
                }
                else if(at.pib() == 1) {

                  // If necessary, calculate drift force
                  if(vp.lambda()<1.0 && vp.lambda()>0.0) {
                    real xi = 0.0;
                    for (std::vector<Particle*>::iterator it5 = atList.begin();
                         it5 != atList.end(); ++it5) {
                      Particle &at2 = **it5;
                      if(at2.pib() != 1) {
                        xi += at2.modepos().sqr()*Eigenvalues[at2.pib()-1];
                      }
                    }
                    std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                    Real3D pa = **it2;
                    Real3D mindriftforce(0.0, 0.0, 0.0);
                    verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                    real min1sq = 0.0;
                    if(verletList->getAdrRegionType()) {
                      min1sq = mindriftforce.sqr();
                    }
                    else {
                      min1sq = mindriftforce[0]*mindriftforce[0];
                    }
                    ++it2;
                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                      pa = **it2;
                      Real3D driftforce(0.0, 0.0, 0.0);
                      verletList->getSystem()->bc->getMinimumImageVector(driftforce, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                      real distsq1 = 0.0;
                      if(verletList->getAdrRegionType()) {
                        distsq1 = driftforce.sqr();
                      }
                      else {
                        distsq1 = driftforce[0]*driftforce[0];
                      }
                      if (distsq1 < min1sq) {
                        min1sq = distsq1;
                        mindriftforce = driftforce;
                      }
                    }
                    min1sq = sqrt(min1sq);
                    if(verletList->getAdrRegionType()) {
                      mindriftforce = (1.0/min1sq)*mindriftforce;
                      mindriftforce *= 0.5*xi*(clmassmultiplier-1.0)*(1.0/sqrt(ntrotter))*vp.lambdaDeriv()*vp.mass()*omega2;
                      at.forcem() += mindriftforce;
                    }
                    else {
                      real mindriftforceX = (1.0/min1sq)*mindriftforce[0];
                      mindriftforceX *= 0.5*xi*(clmassmultiplier-1.0)*(1.0/sqrt(ntrotter))*vp.lambdaDeriv()*vp.mass()*omega2;
                      Real3D driftforceadd(mindriftforceX,0.0,0.0);
                      at.forcem() += driftforceadd;
                    }
                  }

                }
                else {
                  throw std::runtime_error("at.pib() outside of trotter range in calcForcesF function (PIAdressIntegrator)!");
                }
              }
            }
            else {
              std::stringstream ss;
              ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
              throw std::runtime_error(ss.str());
            }
          }

        }
        else {

          if(it2 != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it2->second;
            for (std::vector<Particle*>::iterator it3 = atList.begin();
                 it3 != atList.end(); ++it3) {
              Particle &at = **it3;
              if(at.pib() > 1 && at.pib() <= ntrotter) {
                at.forcem() -= at.modepos()*omega2*vp.varmass()*Eigenvalues[at.pib()-1];
              }
              else if(at.pib() == 1) {

                // If necessary, calculate drift force
                if(vp.lambda()<1.0 && vp.lambda()>0.0) {
                  real xi = 0.0;
                  for (std::vector<Particle*>::iterator it5 = atList.begin();
                       it5 != atList.end(); ++it5) {
                    Particle &at2 = **it5;
                    if(at2.pib() != 1) {
                      xi += at2.modepos().sqr()*Eigenvalues[at2.pib()-1];
                    }
                  }
                  std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                  Real3D pa = **it2;
                  Real3D mindriftforce(0.0, 0.0, 0.0);
                  verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                  real min1sq = 0.0;
                  if(verletList->getAdrRegionType()) {
                    min1sq = mindriftforce.sqr();
                  }
                  else {
                    min1sq = mindriftforce[0]*mindriftforce[0];
                  }
                  ++it2;
                  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                    pa = **it2;
                    Real3D driftforce(0.0, 0.0, 0.0);
                    verletList->getSystem()->bc->getMinimumImageVector(driftforce, (1.0/sqrt(ntrotter))*at.modepos(), pa);
                    real distsq1 = 0.0;
                    if(verletList->getAdrRegionType()) {
                      distsq1 = driftforce.sqr();
                    }
                    else {
                      distsq1 = driftforce[0]*driftforce[0];
                    }
                    if (distsq1 < min1sq) {
                      min1sq = distsq1;
                      mindriftforce = driftforce;
                    }
                  }
                  min1sq = sqrt(min1sq);
                  if(verletList->getAdrRegionType()) {
                    mindriftforce = (1.0/min1sq)*mindriftforce;
                    mindriftforce *= 0.5*xi*(clmassmultiplier-1.0)*(1.0/sqrt(ntrotter))*vp.lambdaDeriv()*vp.mass()*omega2;
                    at.forcem() += mindriftforce;
                  }
                  else {
                    real mindriftforceX = (1.0/min1sq)*mindriftforce[0];
                    mindriftforceX *= 0.5*xi*(clmassmultiplier-1.0)*(1.0/sqrt(ntrotter))*vp.lambdaDeriv()*vp.mass()*omega2;
                    Real3D driftforceadd(mindriftforceX,0.0,0.0);
                    at.forcem() += driftforceadd;
                  }
                }

              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in calcForcesF function (PIAdressIntegrator)!");
              }
            }
          }
          else {
            std::stringstream ss;
            ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
            throw std::runtime_error(ss.str());
          }

        }
      }
    }


    void PIAdressIntegrator::updateForces(int f) {
      initForces(f);

      storage::Storage& storage = *getSystemRef().storage;
      if (f==2 || f==3) {
        storage.updateGhosts();
        if (KTI == false) {
          setWeights();
        }
      }

      if (f==1) {
        calcForcesF();
      }
      else if (f==2) {
        calcForcesM();
      }
      else if (f==3) {
        calcForcesS();
      }
      else {
        throw std::runtime_error("updateForces routine in PIAdressIntegrator received wrong integer.");
      }

      if (f==2 || f==3) { // Not necessary for internal ring forces
        storage.collectGhostForces();
        distributeForces();
        transForces();
      }
    }


    void PIAdressIntegrator::distributeForces() { // distribute forces from atomistic particles (path integral beads) to CG particles (physical atoms)
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);

        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;
          Real3D vpfm = vp.force();
          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            at.force() += (1.0/ntrotter)*vpfm;
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }

      }
    }


    void PIAdressIntegrator::initForces(int f) {
      System& system = getSystemRef();

      if (f==1) { // If mode forces, we do not care about ghosts, CG particles, and real forces
        // real atomistic particles
        CellList localCells = system.storage->getLocalCells();
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
             it != adrATparticles.end(); ++it) {
          it->forcem() = 0.0;
        }
      }

      else if (f==2 || f==3) {
        // all CG particles
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
          cit->force() = 0.0;
          cit->forcem() = 0.0;
        }
        // real atomistic particles
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
             it != adrATparticles.end(); ++it) {
          it->force() = 0.0;
          it->forcem() = 0.0;
        }
        // atomistic ghost particles
        typedef std::list<ParticleList> ParticleListAdr;
        ParticleListAdr& adrATparticlesG = system.storage->getAdrATParticlesG();
        for (ParticleListAdr::iterator it = adrATparticlesG.begin();
             it != adrATparticlesG.end(); ++it) {
          for (ParticleList::iterator it2 = it->begin();
               it2 != it->end(); ++it2) {
            it2->force() = 0.0;
            it2->forcem() = 0.0;
          }
        }
      }

      else {
        throw std::runtime_error("initForces routine in PIAdressIntegrator received wrong integer.");
      }
    }


    real PIAdressIntegrator::computeKineticEnergy() {
      real esum = 0.0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;
          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            if(at.pib() == 1) {
              esum += at.modemom().sqr() * ntrotter / vp.mass();
            }
            else if(at.pib() > 1 && at.pib() <= ntrotter) {
              if((speedup == true) && (vp.lambda() < 0.000000001)) {
                continue;
              } else {
                if(constkinmass == false) {
                  if(realkinmass == false) {
                    esum += at.modemom().sqr() * ntrotter / (vp.varmass()*CMDparameter*Eigenvalues[at.pib()-1]);
                  }
                  else {
                    esum += at.modemom().sqr() * ntrotter / (vp.varmass()*CMDparameter);
                  }
                } else {
                  if(realkinmass == false) {
                    esum += at.modemom().sqr() * ntrotter / (vp.mass()*CMDparameter*Eigenvalues[at.pib()-1]);
                  }
                  else {
                    esum += at.modemom().sqr() * ntrotter / (vp.mass()*CMDparameter);
                  }
                }
              }
            }
            else {
              throw std::runtime_error("at.pib() outside of trotter range in computeKineticEnergy routine (PIAdressIntegrator).");
            }
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
      esum *= 0.5;
      real esumtotal;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());
      return esumtotal;
    }


    real PIAdressIntegrator::computeRingEnergy() { // internal ring energy based on mode positions
      real esum = 0.0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;
          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            if(at.pib() == 1) {
              continue;
            }
            else if(at.pib() > 1 && at.pib() <= ntrotter) {
              esum += at.modepos().sqr() * omega2 *  (vp.varmass()*Eigenvalues[at.pib()-1]);
            }
            else {
              throw std::runtime_error("at.pib() outside of trotter range in computeRingEnergy routine (PIAdressIntegrator).");
            }
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
      esum *= 0.5;
      real esumtotal;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());
      return esumtotal;
    }


    real PIAdressIntegrator::computeRingEnergyRaw() { // internal ring energy based on real positions
      real esum = 0.0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Particle*> atList;
          atList = it3->second;
          Real3D pos1(0.0,0.0,0.0);
          Real3D posN(0.0,0.0,0.0);
          for (std::vector<Particle*>::iterator it2 = atList.begin();
               it2 != atList.end(); ++it2) {
            Particle &at = **it2;
            if(at.pib() == 1) {
              pos1 = at.position();
              posN = at.position();
            }
            else if(at.pib() > 1 && at.pib() <= ntrotter-1) {
              esum += (at.position()-posN).sqr() * omega2 * vp.varmass();
              posN = at.position();
            }
            else if(at.pib() == ntrotter) {
              esum += (at.position()-posN).sqr() * omega2 * vp.varmass();
              esum += (at.position()-pos1).sqr() * omega2 * vp.varmass();
            }
            else {
              throw std::runtime_error("at.pib() outside of trotter range in computeRingEnergyRaw routine (PIAdressIntegrator).");
            }
          }
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }
      esum *= 0.5;
      real esumtotal;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());
      return esumtotal;
    }


    real PIAdressIntegrator::computeMomentumDrift(int parttype) {
      real esum = 0.0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        if(vp.type() == parttype) {
          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);
          if (it3 != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it3->second;
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;
              if(at.pib() == 1) {
                continue;
              }
              else if(at.pib() > 1 && at.pib() <= ntrotter) {
                if(realkinmass == false) {
                  esum += (clmassmultiplier-1.0)*0.5*vp.mass()*ntrotter*at.modemom().sqr()/(vp.varmass()*vp.varmass()*CMDparameter*Eigenvalues[at.pib()-1]);
                }
                else {
                  esum += (clmassmultiplier-1.0)*0.5*vp.mass()*ntrotter*at.modemom().sqr()/(vp.varmass()*vp.varmass()*CMDparameter);
                }
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in computeMomentumDrift routine (PIAdressIntegrator).");
              }
            }
          }
          else {
            std::stringstream ss;
            ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
            throw std::runtime_error(ss.str());
          }
        }
      }
      real esumtotal;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());
      return esumtotal;
    }


    real PIAdressIntegrator::computePositionDrift(int parttype) {
      real esum = 0.0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        if(vp.type() == parttype) {
          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);
          if (it3 != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it3->second;
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                 it2 != atList.end(); ++it2) {
              Particle &at = **it2;
              if(at.pib() == 1) {
                continue;
              }
              else if(at.pib() > 1 && at.pib() <= ntrotter) {
                esum -= (clmassmultiplier-1.0)*0.5*vp.mass()*omega2*at.modepos().sqr()*Eigenvalues[at.pib()-1];
              }
              else {
                throw std::runtime_error("at.pib() outside of trotter range in computePositionDrift routine (PIAdressIntegrator).");
              }
            }
          }
          else {
            std::stringstream ss;
            ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
            throw std::runtime_error(ss.str());
          }
        }
      }
      real esumtotal;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());
      return esumtotal;
    }


    void PIAdressIntegrator::setWeights() {
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {
          std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
          Real3D pa = **it2;
          Real3D d1(0.0, 0.0, 0.0);
          real min1sq;
          verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
          if (verletList->getAdrRegionType()) { // spherical adress region
            min1sq = d1.sqr();
            ++it2;
            for (; it2 != verletList->getAdrPositions().end(); ++it2) {
              pa = **it2;
              verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
              real distsq1 = d1.sqr();
              if (distsq1 < min1sq) min1sq = distsq1;
            }
          }
          else { //slab-type adress region
            min1sq = d1[0]*d1[0];
            ++it2;
            for (; it2 != verletList->getAdrPositions().end(); ++it2) {
              pa = **it2;
              verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
              real distsq1 = d1[0]*d1[0];
              if (distsq1 < min1sq) min1sq = distsq1;
            }
          }
          real w = weight(min1sq);
          vp.lambda() = w;
          real wDeriv = weightderivative(min1sq);
          vp.lambdaDeriv() = wDeriv;
          vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
        }
        else {
          std::stringstream ss;
          ss << "VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples " << " (" << vp.position() << ").";
          throw std::runtime_error(ss.str());
        }
      }

    }


    real PIAdressIntegrator::weight(real distanceSqr) {
      if (dex2 >= distanceSqr) return 1.0;
      else if (dexdhy2 <= distanceSqr) return 0.0;
      else {
        real argument = sqrt(distanceSqr) - dex;
        return pow(cos(pidhy2 * argument),2.0);
      }
    }


    real PIAdressIntegrator::weightderivative(real distanceSqr) {
      if (dex2 >= distanceSqr) return 0.0;
      else if (dexdhy2 <= distanceSqr) return 0.0;
      else {
        real argument = sqrt(distanceSqr) - dex;
        return -1.0 * pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument);
      }
    }


    void PIAdressIntegrator::setTimeStep(real _dt) {
      dt = _dt;
      dt2 = dt*sStep;
      dt3 = dt*mStep*sStep;
    }


    void PIAdressIntegrator::setmStep(int _mStep) {
      mStep = _mStep;
      dt2 = dt*sStep;
      dt3 = dt*mStep*sStep;
    }


    void PIAdressIntegrator::setsStep(int _sStep) {
      sStep = _sStep;
      dt2 = dt*sStep;
      dt3 = dt*mStep*sStep;
    }


    void PIAdressIntegrator::setNtrotter(int _ntrotter) {
      ntrotter = _ntrotter;
    }


    void PIAdressIntegrator::setSpeedup(bool _speedup) {
      speedup = _speedup;
    }


    void PIAdressIntegrator::setKTI(bool _KTI) {
      KTI = _KTI;
    }


    void PIAdressIntegrator::setPILE(bool _PILE) {
      PILE = _PILE;
    }


    void PIAdressIntegrator::setRealKinMass(bool _realkinmass) {
      realkinmass = _realkinmass;
    }


    void PIAdressIntegrator::setCentroidThermostat(bool _centroidthermostat) {
      centroidthermostat = _centroidthermostat;
    }


    void PIAdressIntegrator::setConstKinMass(bool _constkinmass) {
      constkinmass = _constkinmass;
    }


    void PIAdressIntegrator::setTemperature(real _temperature) {
      temperature = _temperature;
      setOmegaSquared();
    }


    void PIAdressIntegrator::setOmegaSquared() {
      real hbar = 1.054571726;
      real Nav = 6.02214129;
      real kb = 1.3806488;
      omega2 = ntrotter*temperature*temperature*kb*kb*Nav*1.66053892/(hbar*hbar*0.00831451*0.00831451*1000.0);
    }


    void PIAdressIntegrator::setGamma(real _gamma) {
      gamma = _gamma;
    }


    void PIAdressIntegrator::setCMDparameter(real _CMDparameter) {
      CMDparameter = _CMDparameter;
    }


    void PIAdressIntegrator::setPILElambda(real _PILElambda) {
      PILElambda = _PILElambda;
    }


    void PIAdressIntegrator::setClmassmultiplier(real _clmassmultiplier) {
      clmassmultiplier = _clmassmultiplier;
    }


    void PIAdressIntegrator::setVerletList(shared_ptr<VerletListAdress> _verletList) {
      if (!_verletList) {
        throw std::invalid_argument("No Verletlist given in PIAdressIntegrator::setVerletList.");
      }
      verletList = _verletList;
    }


    void PIAdressIntegrator::transp() { // We need both the regular and the transposed normal mode transformation matrix
      if (Tvectors.size() != ntrotter)
      {
        throw std::runtime_error("Number of Eigenvectors not equal to number of Trotter beads!");
      }
      std::vector<real> tmpvec;
      for (size_t i = 0; i < ntrotter; ++i)
      {
        for (size_t j = 0; j < ntrotter; ++j)
        {
          tmpvec.push_back( Tvectors.at(j).at(i) );
        }
        Eigenvectors.push_back(tmpvec);
        tmpvec.clear();
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void PIAdressIntegrator::registerPython() {

      using namespace espressopp::python;

      class_<PIAdressIntegrator, bases<MDIntegrator>, boost::noncopyable >
      ("integrator_PIAdressIntegrator", init< shared_ptr<System>, shared_ptr<VerletListAdress> >())
      .def("setTimeStep", &PIAdressIntegrator::setTimeStep)
      .def("getTimeStep", &PIAdressIntegrator::getTimeStep)
      .def("setmStep", &PIAdressIntegrator::setmStep)
      .def("getmStep", &PIAdressIntegrator::getmStep)
      .def("setsStep", &PIAdressIntegrator::setsStep)
      .def("getsStep", &PIAdressIntegrator::getsStep)
      .def("setNtrotter", &PIAdressIntegrator::setNtrotter)
      .def("getNtrotter", &PIAdressIntegrator::getNtrotter)
      .def("setTemperature", &PIAdressIntegrator::setTemperature)
      .def("getTemperature", &PIAdressIntegrator::getTemperature)
      .def("setGamma", &PIAdressIntegrator::setGamma)
      .def("getGamma", &PIAdressIntegrator::getGamma)
      .def("setCMDparameter", &PIAdressIntegrator::setCMDparameter)
      .def("getCMDparameter", &PIAdressIntegrator::getCMDparameter)
      .def("setPILElambda", &PIAdressIntegrator::setPILElambda)
      .def("getPILElambda", &PIAdressIntegrator::getPILElambda)
      .def("setClmassmultiplier", &PIAdressIntegrator::setClmassmultiplier)
      .def("getClmassmultiplier", &PIAdressIntegrator::getClmassmultiplier)
      .def("setSpeedup", &PIAdressIntegrator::setSpeedup)
      .def("getSpeedup", &PIAdressIntegrator::getSpeedup)
      .def("setKTI", &PIAdressIntegrator::setKTI)
      .def("getKTI", &PIAdressIntegrator::getKTI)
      .def("setCentroidThermostat", &PIAdressIntegrator::setCentroidThermostat)
      .def("getCentroidThermostat", &PIAdressIntegrator::getCentroidThermostat)
      .def("setPILE", &PIAdressIntegrator::setPILE)
      .def("getPILE", &PIAdressIntegrator::getPILE)
      .def("setRealKinMass", &PIAdressIntegrator::setRealKinMass)
      .def("getRealKinMass", &PIAdressIntegrator::getRealKinMass)
      .def("setConstKinMass", &PIAdressIntegrator::setConstKinMass)
      .def("getConstKinMass", &PIAdressIntegrator::getConstKinMass)
      .def("setVerletList", &PIAdressIntegrator::setVerletList)
      .def("getVerletList", &PIAdressIntegrator::getVerletList)
      .def("addEVcomponent", &PIAdressIntegrator::addEVcomponent)
      .def("addTransposedEigenVector", &PIAdressIntegrator::addTransposedEigenVector)
      .def("addEigenValues", &PIAdressIntegrator::addEigenValues)
      .def("transp", &PIAdressIntegrator::transp)
      .def("computeKineticEnergy", &PIAdressIntegrator::computeKineticEnergy)
      .def("computeRingEnergy", &PIAdressIntegrator::computeRingEnergy)
      .def("computeRingEnergyRaw", &PIAdressIntegrator::computeRingEnergyRaw)
      .def("computeMomentumDrift", &PIAdressIntegrator::computeMomentumDrift)
      .def("computePositionDrift", &PIAdressIntegrator::computePositionDrift)
      .def("getVerletlistBuilds", &PIAdressIntegrator::getVerletlistBuilds)
      .add_property("timestep", &PIAdressIntegrator::getTimeStep, &PIAdressIntegrator::setTimeStep)
      .add_property("mSteps", &PIAdressIntegrator::getmStep, &PIAdressIntegrator::setmStep)
      .add_property("sSteps", &PIAdressIntegrator::getsStep, &PIAdressIntegrator::setsStep)
      .add_property("nTrotter", &PIAdressIntegrator::getNtrotter, &PIAdressIntegrator::setNtrotter)
      .add_property("temperature", &PIAdressIntegrator::getTemperature, &PIAdressIntegrator::setTemperature)
      .add_property("gamma", &PIAdressIntegrator::getGamma, &PIAdressIntegrator::setGamma)
      .add_property("CMDparameter", &PIAdressIntegrator::getCMDparameter, &PIAdressIntegrator::setCMDparameter)
      .add_property("CLmassmultiplier", &PIAdressIntegrator::getClmassmultiplier, &PIAdressIntegrator::setClmassmultiplier)
      .add_property("PILElambda", &PIAdressIntegrator::getPILElambda, &PIAdressIntegrator::setPILElambda)
      .add_property("speedup", &PIAdressIntegrator::getSpeedup, &PIAdressIntegrator::setSpeedup)
      .add_property("KTI", &PIAdressIntegrator::getKTI, &PIAdressIntegrator::setKTI)
      .add_property("centroidThermostat", &PIAdressIntegrator::getCentroidThermostat, &PIAdressIntegrator::setCentroidThermostat)
      .add_property("PILE", &PIAdressIntegrator::getPILE, &PIAdressIntegrator::setPILE)
      .add_property("realKinMass", &PIAdressIntegrator::getRealKinMass, &PIAdressIntegrator::setRealKinMass)
      .add_property("constKinMass", &PIAdressIntegrator::getConstKinMass, &PIAdressIntegrator::setConstKinMass)
      .add_property("verletList", &PIAdressIntegrator::getVerletList, &PIAdressIntegrator::setVerletList)
      .add_property("verletlistBuilds", &PIAdressIntegrator::getVerletlistBuilds)
      ;
    }
  }
}
