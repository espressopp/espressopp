#include <iomanip>
#include "python.hpp"
#include "VelocityVerlet.hpp"

/* TODO should be removed after signals will be tested
#include "Langevin.hpp"
#include "BerendsenBarostat.hpp"
#include "BerendsenThermostat.hpp"
#include "LangevinBarostat.hpp"
#include "Isokinetic.hpp"
#include "StochasticVelocityRescaling.hpp"
#include "FixPositions.hpp"
*/

#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espresso {

  using namespace std;
  
  namespace integrator {

    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    VelocityVerlet::VelocityVerlet(shared_ptr< System > system) : MDIntegrator(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");

      resortFlag = true;
      maxDist    = 0.0;
    }

    /*****************************************************************************/

    VelocityVerlet::~VelocityVerlet()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerlet");
    }

    // TODO should be removed after signals will be tested
    /*****************************************************************************
    void VelocityVerlet::setLangevin(shared_ptr< Langevin > _langevin)
    {
      LOG4ESPP_INFO(theLogger, "set Langevin thermostat");
      langevin = _langevin;
    }

    void VelocityVerlet::setIsokinetic(shared_ptr< Isokinetic > _isokinetic)
    {
      LOG4ESPP_INFO(theLogger, "set Isokinetic thermostat");
      isokinetic = _isokinetic;
    }

    void VelocityVerlet::setStochasticVelocityRescaling(shared_ptr< StochasticVelocityRescaling > _stochasticVelocityRescaling)
    {
      LOG4ESPP_INFO(theLogger, "set StochasticVelocityRescaling thermostat");
      stochasticVelocityRescaling = _stochasticVelocityRescaling;
    }

    void VelocityVerlet::setBerendsenBarostat(shared_ptr< BerendsenBarostat > _berendsenP){
      LOG4ESPP_INFO(theLogger, "set Berendsen barostat");
      berendsenBarostat = _berendsenP;
    }

    void VelocityVerlet::setBerendsenThermostat(shared_ptr< BerendsenThermostat > _berendsenT){
      LOG4ESPP_INFO(theLogger, "set Berendsen thermostat");
      berendsenThermostat = _berendsenT;
    }

    void VelocityVerlet::setLangevinBarostat(shared_ptr< LangevinBarostat > _langevinBarostat){
      LOG4ESPP_INFO(theLogger, "set Hoover style barostat. Langevin dynamics");
      langevinBarostat = _langevinBarostat;
    }
    
    void VelocityVerlet::setFixPositions(shared_ptr <class FixPositions> _fixPositions) {
    	LOG4ESPP_INFO(theLogger, "set fixPositions to fix positions of a group of particles");
    	fixPositions = _fixPositions;
    }
    *****************************************************************************/

    void VelocityVerlet::run(int nsteps)
    {
      VT_TRACER("run");

      int nResorts = 0;

      real time;

      timeIntegrate.reset();

      resetTimers();

      System& system = getSystemRef();

      storage::Storage& storage = *system.storage;


      // signal
      runInit();

      /* -- use signal instead
      if (langevin) langevin->initialize(dt);

      if (langevinBarostat) langevinBarostat->initialize(dt);
      
      if (berendsenBarostat) berendsenBarostat->initialize(dt);
      if (berendsenThermostat) berendsenThermostat->initialize(dt);
      */

      
      // no more needed: setUp();

      // Before start make sure that particles are on the right processor

      if (resortFlag) {
        VT_TRACER("resort");
        // time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        LOG4ESPP_INFO(theLogger, "particles resort");
        maxDist = 0.0;
        resortFlag = false;
        // timeResort += timeIntegrate.getElapsedTime();
      }

      bool recalcForces = true;  // TODO: more intelligent

      if (recalcForces) {

        LOG4ESPP_INFO(theLogger, "recalc Forces");

        // signal
        recalc1();

        //if (langevin) langevin->heatUp();
        //if (langevinBarostat) langevinBarostat->heatUp();   do it needs heatUp analog?
        
        updateForces();
        if (LOG4ESPP_DEBUG_ON(theLogger)) {
            // printForces(false);   // forces are reduced to real particles
        }

        // signal
        recalc2();

        //if (langevin) langevin->coolDown();
        //if (langevinBarostat) langevinBarostat->coolDown();  the same question like for heatUp
      }

      LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
      real skinHalf = 0.5 * system.skin;

      for (int i = 0; i < nsteps; i++) {

        LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        //saveOldPos(); // save particle positions needed for constraints

        // signal
        befIntP();


        time = timeIntegrate.getElapsedTime();
        maxDist += integrate1();
        timeInt1 += timeIntegrate.getElapsedTime() - time;


        // signal
        aftIntP();


        LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag) {
            VT_TRACER("resort1");
            time = timeIntegrate.getElapsedTime();
            LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
            storage.decompose();
            maxDist  = 0.0;
            resortFlag = false;
            nResorts ++;
            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        updateForces();


        // signal
        befIntV();


        time = timeIntegrate.getElapsedTime();
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;
        

        // signal
        aftIntV();


        /* -- use signals instead
        if (berendsenThermostat) berendsenThermostat->thermostat(); // adjust temperature
        if (berendsenBarostat) berendsenBarostat->barostat(); // adjust pressure
        if (isokinetic) isokinetic->rescaleVelocities(); // scale all particle velocities to match isokinetic temperature
        if (stochasticVelocityRescaling) stochasticVelocityRescaling->rescaleVelocities(); // scale all particle velocities according to a canonical distribution
        */

      }

      timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);

      LOG4ESPP_INFO(theLogger, "finished run");

      // ToDo: print Timers only if INFO is enabled
      //printTimers();
    }

    void VelocityVerlet::resetTimers() {
      timeForce  = 0.0;
      for(int i = 0; i < 100; i++)
        timeForceComp[i] = 0.0;
      timeComm1  = 0.0;
      timeComm2  = 0.0;
      timeInt1   = 0.0;
      timeInt2   = 0.0;
      timeResort = 0.0;
    }

    using namespace boost::python;

    static object wrapGetTimers(class VelocityVerlet* obj) {
      real tms[10];
      obj->loadTimers(tms);
      return make_tuple(tms[0],
                        tms[1],
                        tms[2],
                        tms[3],
                        tms[4],
                        tms[5],
                        tms[6],
                        tms[7],
                        tms[8],
                        tms[9]);
    }

    void VelocityVerlet::loadTimers(real t[10]) {
      t[0] = timeRun;
      t[1] = timeForceComp[0];
      t[2] = timeForceComp[1];
      t[3] = timeForceComp[2];
      t[4] = timeComm1;
      t[5] = timeComm2;
      t[6] = timeInt1;
      t[7] = timeInt2;
      t[8] = timeResort;
      t[9] = timeLost;
    }

    void VelocityVerlet::printTimers() {

      using namespace std;

      real pct;

      cout << endl;
      cout << "run = " << setiosflags(ios::fixed) << setprecision(1) << timeRun << endl;
      pct = 100.0 * (timeForceComp[0] / timeRun);
      cout << "pair (%) = " << timeForceComp[0] << " (" << pct << ")" << endl;
      pct = 100.0 * (timeForceComp[1] / timeRun);
      cout << "FENE (%) = " << timeForceComp[1] << " (" << pct << ")" << endl;
      pct = 100.0 * (timeForceComp[2] / timeRun);
      cout << "angle (%) = " << timeForceComp[2] << " (" << pct << ")" << endl;
      pct = 100.0 * (timeComm1 / timeRun);
      cout << "comm1 (%) = " << timeComm1 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeComm2 / timeRun);
      cout << "comm2 (%) = " << timeComm2 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeInt1 / timeRun);
      cout << "int1 (%) = " << timeInt1 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeInt2 / timeRun);
      cout << "int2 (%) = " << timeInt2 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeResort / timeRun);
      cout << "resort (%) = " << timeResort << " (" << pct << ")" << endl;
      pct = 100.0 * (timeLost / timeRun);
      cout << "other (%) = " << timeLost << " (" << pct << ")" << endl;
      cout << endl;
    }

    /*****************************************************************************/

    real VelocityVerlet::integrate1()
    {
      System& system = getSystemRef();
      real half_dt = 0.5 * dt; 


      /* -- use signal instead
      if (langevinBarostat){
        // update the volume V(t+0.5*dt)=V(t)+dt/2. * V'
        langevinBarostat->updVolume( dt );
        // update the local barostat momentum pe(t+0.5*dt)=pe(t)+dt/2. * pe'
        langevinBarostat->updVolumeMomentum( half_dt );
      }
      */
      
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells

      int count = 0;

      real maxSqDist = 0.0; // maximal square distance a particle moves

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

        real sqDist = 0.0;

        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << 
                ", pos = " << cit->position() <<
                ", v = " << cit->velocity() << 
                ", f = " << cit->force());

        /* more precise for DEBUG:

        printf("Particle %d, pos = %16.12f %16.12f %16.12f, v = %16.12f, %16.12f %16.12f, f = %16.12f %16.12f %16.12f\n",
            cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2],
                cit->m.v[0], cit->m.v[1], cit->m.v[2],
            cit->f.f[0], cit->f.f[1], cit->f.f[2]);

        */

        real dtfm = half_dt / cit->mass();

        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) 
        cit->velocity() += dtfm * cit->force();

        // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt) 
        Real3D deltaP = cit->velocity();

        
        // use signal instead ?
        /*
        if(langevinBarostat) deltaP += (langevinBarostat->updDisplacement()) * cit->position(); // updDisplacement is just coefficient

        deltaP *= dt;

        if (fixPositions) fixPositions->apply(cit->id(), cit->velocity(), deltaP);

        */

        deltaP *= dt;
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;


        count++;

        maxSqDist = std::max(maxSqDist, sqDist);
      }
      


      // signal
      inIntP(maxSqDist);

      real maxAllSqDist;

      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));
      
      return sqrt(maxAllSqDist);
    }

    /*****************************************************************************/

    void VelocityVerlet::integrate2()
    {
      System& system = getSystemRef();

      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells

      real half_dt = 0.5 * dt; 
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

        real dtfm = half_dt / cit->mass();

        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */

        cit->velocity() += dtfm * cit->force();
      }
      

      /* -- use signal instead
      if (langevinBarostat){
        // update the local momentum pe(t+dt)=pe(t+0.5*dt)+dt/2. * pe'
        langevinBarostat->updVolumeMomentum( half_dt );
        // update the volume V(t+dt)=V(t+0.5*dt)+dt/2. * V'
        langevinBarostat->updVolume( dt );
      }*/

      
      step++;
    }

    /*************DO NOT NEED, PROBABLY SHOULD BE DELETED*************************

    void VelocityVerlet::setUp()
    {
      System& system = getSystemRef();

      const InteractionList& srIL = system.shortRangeInteractions;

      maxCut = 0.0;

      for (size_t j = 0; j < srIL.size(); j++) {
        real cut = srIL[j]->getMaxCutoff();
        maxCut = std::max(maxCut, cut);
      }

      LOG4ESPP_INFO(theLogger, "maximal cutoff = " << maxCut);
    }

    /*****************************************************************************/

    void VelocityVerlet::calcForces()
    {
      VT_TRACER("forces");

      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

      // signal
      aftInitF();


      System& sys = getSystemRef();

      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {

	LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i 
                                  << " of " << srIL.size());

        real time;
        time = timeIntegrate.getElapsedTime();
        srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
      }
    }

    void VelocityVerlet::updateForces()
    {
      real time;
      storage::Storage& storage = *getSystemRef().storage;

      time = timeIntegrate.getElapsedTime();
      { 
        VT_TRACER("commF");
        storage.updateGhosts();
      }
      timeComm1 += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      calcForces();
      timeForce += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      {
        VT_TRACER("commR");
        storage.collectGhostForces();
      }
      timeComm2 += timeIntegrate.getElapsedTime() - time;

      // signal
      aftCalcF();

      //if (langevin) langevin->thermalize();
      //if (langevinBarostat) langevinBarostat->updForces();
    }

    /*****************************************************************************/

    void VelocityVerlet::initForces()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();

      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
      }
    }

    /*****************************************************************************/

    /*
    System& getSystem() {
        return getsystemRef();
    }*/


    void VelocityVerlet::printForces(bool withGhosts)
    {
      // print forces of real + ghost particles

      System& system = getSystemRef();

      CellList cells;

      if (withGhosts) {
	cells = system.storage->getLocalCells();
	LOG4ESPP_DEBUG(theLogger, "local forces");
      } else {
	cells = system.storage->getRealCells();
	LOG4ESPP_DEBUG(theLogger, "real forces");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

	LOG4ESPP_DEBUG(theLogger, 
		       "Particle " << cit->id()
		       << ", force = " << cit->force());
      }
    }

    /*****************************************************************************/

    void VelocityVerlet::printPositions(bool withGhosts)
    {
      // print positions of real + ghost particles

      System& system = getSystemRef();

      CellList cells;

      if (withGhosts) {
	cells = system.storage->getLocalCells();
	LOG4ESPP_DEBUG(theLogger, "local positions");
      } else {
	cells = system.storage->getRealCells();
	LOG4ESPP_DEBUG(theLogger, "real positions");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

	LOG4ESPP_DEBUG(theLogger, 
		       "Particle " << cit->id()
		       << ", position = " << cit->position());
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerlet::registerPython() {

      using namespace espresso::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerlet, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerlet", init< shared_ptr<System> >())
        /* TODO should be removed after signals will be tested
        .add_property("langevin", &VelocityVerlet::getLangevin, &VelocityVerlet::setLangevin)
        .add_property("berendsenBarostat", &VelocityVerlet::getBerendsenBarostat, &VelocityVerlet::setBerendsenBarostat)
        .add_property("berendsenThermostat", &VelocityVerlet::getBerendsenThermostat, &VelocityVerlet::setBerendsenThermostat)
        .add_property("langevinBarostat", &VelocityVerlet::getLangevinBarostat, &VelocityVerlet::setLangevinBarostat)
        .add_property("isokinetic", &VelocityVerlet::getIsokinetic, &VelocityVerlet::setIsokinetic)
        .add_property("stochasticVelocityRescaling", &VelocityVerlet::getStochasticVelocityRescaling, &VelocityVerlet::setStochasticVelocityRescaling)
        .add_property("fixpositions", &VelocityVerlet::getFixPositions, &VelocityVerlet::setFixPositions)
         */
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerlet::resetTimers)
        ;
    }
  }
}
