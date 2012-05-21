// ESPP_CLASS
#ifndef _INTEGRATOR_VELOCITYVERLET_HPP
#define _INTEGRATOR_VELOCITYVERLET_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include "esutil/Timer.hpp"
#include <boost/signals2.hpp>

namespace espresso {
  namespace integrator {

    /** Velocity Verlet Integrator */
    class VelocityVerlet : public MDIntegrator {

      public:

        VelocityVerlet(shared_ptr<class espresso::System> system);

        ~VelocityVerlet();

        // set & get termostat (Langevin)
        void setLangevin(shared_ptr<class Langevin> langevin);
        shared_ptr<class Langevin> getLangevin() { return langevin; }

        void setStochasticVelocityRescaling(shared_ptr<class StochasticVelocityRescaling> stochasticVelocityRescaling);
        shared_ptr<class StochasticVelocityRescaling> getStochasticVelocityRescaling() { return stochasticVelocityRescaling; }

        void setIsokinetic(shared_ptr<class Isokinetic> isokinetic);
        shared_ptr<class Isokinetic> getIsokinetic() { return isokinetic; }
        
        // set & get barostat (Berendsen)
        void setBerendsenBarostat(shared_ptr<class BerendsenBarostat> berendsenBarostat);
        shared_ptr<class BerendsenBarostat> getBerendsenBarostat() { return berendsenBarostat; }
        
        // set & get thermostat (Berendsen)
        void setBerendsenThermostat(shared_ptr<class BerendsenThermostat> berendsenThermostat);
        shared_ptr<class BerendsenThermostat> getBerendsenThermostat() { return berendsenThermostat; }
        
        // set & get barostat (Langevin-Hoover)
        void setLangevinBarostat(shared_ptr<class LangevinBarostat> langevinBarostat);
        shared_ptr<class LangevinBarostat> getLangevinBarostat() { return langevinBarostat; }
        
        // set & get FixPositions
        void setFixPositions(shared_ptr <class FixPositions> _fixPositions);
        shared_ptr<class FixPositions> getFixPositions () {return fixPositions; }

        void run(int nsteps);
        
        /** Load timings in array to export to Python as a tuple. */
        void loadTimers(real t[10]);

        void resetTimers();


        // signal used for constraints
        //boost::signals2::signal0 <void> saveOldPos;
		//boost::signals2::signal0 <void> applyPosConst;
		//boost::signals2::signal0 <void> applyVelConst;


        /* -- moved to superclass
        // signals to extend the integrator
        boost::signals2::signal0 <void> runInit; // initialization of run()
        boost::signals2::signal0 <void> recalc1; // inside recalc, before updateForces()
        boost::signals2::signal0 <void> recalc2; // inside recalc, after  updateForces()
        boost::signals2::signal0 <void> befIntP; // before integrate1()
        boost::signals2::signal1 <void, real&> inIntP; // inside end of integrate1()
        boost::signals2::signal0 <void> aftIntP; // after  integrate1()
        boost::signals2::signal0 <void> aftInitF; // after initForces()
        boost::signals2::signal0 <void> aftCalcF; // after calcForces()
        boost::signals2::signal0 <void> befIntV; // before integrate2()
        boost::signals2::signal0 <void> aftIntV; // after  integrate2()
        */

        //System& getSystem();


        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        bool resortFlag;  //!< true implies need for resort of particles
        real maxDist;

        real maxCut;


        shared_ptr< class Langevin > langevin;  //!< Langevin thermostat if available

        shared_ptr< class Isokinetic > isokinetic;  //!< Isokinetic thermostat if available
        shared_ptr< class StochasticVelocityRescaling > stochasticVelocityRescaling;  //!< Stochastic velocity rescaling thermostat if available

        /** Method updates particle positions and velocities.
            \return maximal square distance a particle has moved.
        */

        shared_ptr< class BerendsenBarostat > berendsenBarostat;  //!< Berendsen barostat if available
        shared_ptr< class BerendsenThermostat> berendsenThermostat;  //!< Berendsen thermostat if available
        
        shared_ptr< class LangevinBarostat > langevinBarostat;  //!< Langevin-Hoover barostat if available
        
        shared_ptr< class FixPositions > fixPositions; // fix positions of a group of particles

        real integrate1();

        void integrate2();

        void initForces();

        void updateForces();

        void calcForces();

        void printPositions(bool withGhost);

        void printForces(bool withGhost);

        void setUp();   //!< set up for a new run

        void printTimers();

        esutil::WallTimer timeIntegrate;  //!< used for timing

        // variables that keep time information about different phases
        real timeRun;
        real timeLost;
        real timeForce;
        real timeForceComp[100];
        real timeComm1;
        real timeComm2;
        real timeInt1;
        real timeInt2;
        real timeResort;
    };
  }
}

#endif
