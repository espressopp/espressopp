#include <iomanip>
#include "python.hpp"
#include "VelocityVerletAdress.hpp"

#include "Langevin.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"


namespace espresso {
  namespace integrator {

    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    VelocityVerletAdress::VelocityVerletAdress(shared_ptr< System > system) : MDIntegrator(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerletAdress");

      resortFlag = true;
    }

    /*****************************************************************************/

    VelocityVerletAdress::~VelocityVerletAdress()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerletAdress");
    }

    /*****************************************************************************/

    void VelocityVerletAdress::setLangevin(shared_ptr< Langevin > _langevin)
    {
      LOG4ESPP_INFO(theLogger, "set Langevin thermostat");
      langevin = _langevin;
    }

    /*****************************************************************************/

    void VelocityVerletAdress::run(int nsteps) {
      int nResorts = 0;

      real time;

      timeIntegrate.reset();

      resetTimers();

      System& system = getSystemRef();

      storage::Storage& storage = *system.storage;

      if (langevin) langevin->initialize(dt);

      // no more needed: setUp();

      // Before start make sure that particles are on the right processor

      real maxDist;

      if (resortFlag) {
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
        if (langevin) langevin->heatUp();
        updateForces();
        if (LOG4ESPP_DEBUG_ON(theLogger)) {
            // printForces(false);   // forces are reduced to real particles
        }
        if (langevin) langevin->coolDown();
      }

      LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
      real skinHalf = 0.5 * system.skin;

      for (int i = 0; i < nsteps; i++) {

        LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        saveOldPos(); // save particle positions needed for constraints

        time = timeIntegrate.getElapsedTime();
        maxDist += integrate1();
        timeInt1 += timeIntegrate.getElapsedTime() - time;

        applyConstraints(); // apply constraints

        LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag) {
            time = timeIntegrate.getElapsedTime();
            LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
            //std::cout << " resort particles\n";
            storage.decompose();
            maxDist  = 0.0;
            resortFlag = false;
            nResorts ++;
            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        updateForces();

        time = timeIntegrate.getElapsedTime();
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;
      }

      timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);

      LOG4ESPP_INFO(theLogger, "finished run");

      // ToDo: print Timers only if INFO is enabled
      //printTimers();
    }

    void VelocityVerletAdress::resetTimers() {
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

    static object wrapGetTimers(class VelocityVerletAdress* obj) {
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

    void VelocityVerletAdress::loadTimers(real t[10]) {
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

    void VelocityVerletAdress::printTimers() {

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

    real VelocityVerletAdress::integrate1() {
      System& system = getSystemRef();
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

            real dtfm = 0.5 * dt / cit->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) 
            cit->velocity() += dtfm * cit->force();

            // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt) 
            Real3D deltaP = dt * cit->velocity();
            //std::cout << cit->id() << ": from (" << cit->position() << ")";
            cit->position() += deltaP;
            sqDist += deltaP * deltaP;
            //std::cout << " to (" << cit->position() << ") " << sqrt(sqDist) << "\n";

            count++;

            maxSqDist = std::max(maxSqDist, sqDist);
      }



      // for AdResS
      // propagate real AT particles
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
              it != adrATparticles.end(); it++) {

          real sqDist = 0.0;

          real dtfm = 0.5 * dt / it->mass();

          // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
          it->velocity() += dtfm * it->force();

          // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
          Real3D deltaP = dt * it->velocity();
          //std::cout << "Moving AT " << it->id() << " (" << &(*it) << "): from (" << it->position() << ")";
          it->position() += deltaP;
          sqDist += deltaP * deltaP;
          //std::cout << " to (" << it->position() << ")\n";

          maxSqDist = std::max(maxSqDist, sqDist);
      }


      real maxAllSqDist;

      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, 
                      boost::mpi::maximum<real>());

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));

      return sqrt(maxAllSqDist);
    }

    /*****************************************************************************/

    void VelocityVerletAdress::integrate2()
    {
      System& system = getSystemRef();

      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

        real dtfm = 0.5 * dt / cit->mass();

        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */

        cit->velocity() += dtfm * cit->force();
      }


      // for AdResS
      // propagete real AT particles
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
              it != adrATparticles.end(); ++it) {

          real dtfm = 0.5 * dt / it->mass();

          // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
          it->velocity() += dtfm * it->force();
      }


      step++;
    }

    /*****************************************************************************/

    void VelocityVerletAdress::setUp()
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

    void VelocityVerletAdress::calcForces()
    {
      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

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

    void VelocityVerletAdress::updateForces()
    {
      real time;
      storage::Storage& storage = *getSystemRef().storage;

      time = timeIntegrate.getElapsedTime();
      storage.updateGhosts();
      timeComm1 += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      calcForces();
      timeForce += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      storage.collectGhostForces();
      timeComm2 += timeIntegrate.getElapsedTime() - time;

      if (langevin) langevin->thermalizeAdr();
    }

    /*****************************************************************************/

    void VelocityVerletAdress::initForces()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();

      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
      }


      // for AdResS
      // AT reals
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
            it != adrATparticles.end(); ++it) {
         it->force() = 0.0;

         /*
         if (it->getId() == 12573) {
             std::cout << "resetting force of real  12573 (" << &(*it) << ")\n";
         }*/

      }

      // AT ghosts
      /*
      typedef std::list<Particle> ParticleListAdr;
      ParticleListAdr& adrATparticlesG = system.storage->getAdrATParticlesG();
      for (ParticleListAdr::iterator it = adrATparticlesG.begin();
            it != adrATparticlesG.end(); it++) {
         it->force() = 0.0;
      }*/

      typedef std::list<ParticleList> ParticleListAdr;
      ParticleListAdr& adrATparticlesG = system.storage->getAdrATParticlesG();
      for (ParticleListAdr::iterator it = adrATparticlesG.begin();
              it != adrATparticlesG.end(); ++it) {

          //ParticleList& atgl = *it;

          //std::cout << "atgl size: " << atgl.size() << "\n";

          for (ParticleList::iterator it2 = it->begin();
                  it2 != it->end(); ++it2) {

              it2->force() = 0.0;

              //Particle& atg = *it2;

              /*
              if (atg.getId() == 12573) {
                  std::cout << "resetting force of ghost 12573 (" << &atg << ") pos (" <<
                        atg.getPos() << ")\n";
              }*/

          }

      }


    }

    /*****************************************************************************/

    void VelocityVerletAdress::printForces(bool withGhosts)
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

    void VelocityVerletAdress::printPositions(bool withGhosts)
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

    void VelocityVerletAdress::registerPython() {

      using namespace espresso::python;

      // Note: use noncopyable and no_init for abstract classes
      class_< VelocityVerletAdress, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerletAdress", init< shared_ptr<System> >())
        .add_property("langevin", &VelocityVerletAdress::getLangevin, &VelocityVerletAdress::setLangevin)
        .def("getTimers", &wrapGetTimers)
        ;
    }
  }
}
