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

#ifndef _INTEGRATOR_PIADRESSINTEGRATOR_HPP
#define _INTEGRATOR_PIADRESSINTEGRATOR_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include <boost/signals2.hpp>
#include "VerletListAdress.hpp"

namespace espressopp {
  namespace integrator {

    class PIAdressIntegrator : public MDIntegrator {

      public:
        PIAdressIntegrator(shared_ptr<class espressopp::System> system, shared_ptr<VerletListAdress> _verletList);

        virtual ~PIAdressIntegrator();

        shared_ptr<VerletListAdress> verletList;

        void run(int nsteps);

        void addEVcomponent(real val) { tmpvals.push_back(val); }
        void addTransposedEigenVector() { Tvectors.push_back(tmpvals); tmpvals.clear();}
        void addEigenValues(real val) { Eigenvalues.push_back(val); }
        void transp();

        int getVerletlistBuilds() { return verletlistBuilds; }

        void setTimeStep(real _dt);
        real getTimeStep() { return dt; }

        void setmStep(int _mStep);
        int getmStep() { return mStep; }

        void setsStep(int _sStep);
        int getsStep() { return sStep; }

        void setNtrotter(int _ntrotter);
        int getNtrotter() { return ntrotter; }

        void setTemperature(real _temperature);
        real getTemperature() { return temperature; }

        void setGamma(real _gamma);
        real getGamma() { return gamma; }

        void setCMDparameter(real _CMDparameter);
        real getCMDparameter() { return CMDparameter; }

        void setPILElambda(real _PILElambda);
        real getPILElambda() { return PILElambda; }

        void setClmassmultiplier(real _clmassmultiplier);
        real getClmassmultiplier() { return clmassmultiplier; }

        void setSpeedup(bool _speedup);
        bool getSpeedup() { return speedup; }

        void setKTI(bool _KTI);
        bool getKTI() { return KTI; }

        void setCentroidThermostat(bool _centroidthermostat);
        bool getCentroidThermostat() { return centroidthermostat; }

        void setPILE(bool _PILE);
        bool getPILE() { return PILE; }

        void setRealKinMass(bool _realkinmass);
        bool getRealKinMass() { return realkinmass; }

        void setConstKinMass(bool _constkinmass);
        bool getConstKinMass() { return constkinmass; }

        void setVerletList(shared_ptr<VerletListAdress> _verletList);
        shared_ptr<VerletListAdress> getVerletList() { return verletList; }

        real computeRingEnergy();
        real computeRingEnergyRaw();
        real computeKineticEnergy();

        real computePositionDrift(int parttype);
        real computeMomentumDrift(int parttype);

        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        int sStep;
        int mStep;
        int ntrotter;
        int verletlistBuilds;
        bool resortFlag;
        bool speedup;
        bool KTI;
        bool constkinmass;
        bool realkinmass;
        bool centroidthermostat;
        bool PILE;

        real omega2;
        real clmassmultiplier;

        real CMDparameter;
        real PILElambda;

        real gamma;
        real temperature;

        real maxDist;
        real dt2;
        real dt3;

        real dhy;
        real pidhy2;
        real dex;
        real dex2;
        real dexdhy;
        real dexdhy2;

        std::vector< std::vector<real> > Eigenvectors;
        std::vector< std::vector<real> > Tvectors;
        std::vector< real > Eigenvalues;
        std::vector<real> tmpvals;

        shared_ptr< esutil::RNG > rng;

        void integrateV1(int t, bool doubletime);
        void integrateV2();
        void integrateModePos();
        void OUintegrate();

        void initForces(int f);
        void updateForces(int f);
        void distributeForces();

        void calcForcesF();
        void calcForcesM();
        void calcForcesS();

        void transForces();
        void transPos1();
        void transPos2();
        void transMom1();
        void transMom2();

        void initializeSetup();
        void setWeights();
        void setOmegaSquared();
        real weight(real distanceSqr);
        real weightderivative(real distanceSqr);
    };
  }
}

#endif
