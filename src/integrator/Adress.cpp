#include "python.hpp"
#include "Adress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

  namespace integrator {

    using namespace espresso::iterator;

    Adress::Adress(shared_ptr<System> system)
        :Extension(system){
        LOG4ESPP_INFO(theLogger, "construct Adress");
        type = Extension::Adress;
    }


    Adress::~Adress() {
      LOG4ESPP_INFO(theLogger, "~Adress");
      disconnect();
    }

    void Adress::disconnect(){
        _initForces.disconnect();
        _integrate1.disconnect();
        _integrate2.disconnect();
    }

    void Adress::connect() {

        // connection to after initForces()
        _initForces = integrator->aftInitF.connect(
                boost::bind(&Adress::initForces, this));

        // connection to inside of integrate1()
        _integrate1 = integrator->inIntP.connect(
                boost::bind(&Adress::integrate1, this, _1));

        // connection to after integrate2()
        _integrate2 = integrator->aftIntV.connect(
                boost::bind(&Adress::integrate2, this));
    }


    void Adress::initForces(){

        System& system = getSystemRef();

        // AT reals
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); ++it) {
            it->force() = 0.0;
        }

        // AT ghosts
        typedef std::list<ParticleList> ParticleListAdr;
        ParticleListAdr& adrATparticlesG = system.storage->getAdrATParticlesG();
        for (ParticleListAdr::iterator it = adrATparticlesG.begin();
                it != adrATparticlesG.end(); ++it) {

            for (ParticleList::iterator it2 = it->begin();
                    it2 != it->end(); ++it2) {
                it2->force() = 0.0;
            }
        }
    }



    void Adress::integrate1(real& maxSqDist){

        System& system = getSystemRef();
        real dt = integrator->getTimeStep();

        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); it++) {

            real sqDist = 0.0;
            real dtfm = 0.5 * dt / it->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
            it->velocity() += dtfm * it->force();

            // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
            Real3D deltaP = dt * it->velocity();
            //std::cout << it->id() << ": from (" << it->position() << ")";
            it->position() += deltaP;
            sqDist += deltaP * deltaP;
            //std::cout << " to (" << it->position() << ") " << sqrt(sqDist) << "\n";

            maxSqDist = std::max(maxSqDist, sqDist);
        }

        //std::cout << " " << maxSqDist << "\n";
    }


    void Adress::integrate2() {

        System& system = getSystemRef();
        real dt = integrator->getTimeStep();

        // propagete real AT particles
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); ++it) {

            real dtfm = 0.5 * dt / it->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
            it->velocity() += dtfm * it->force();
        }
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Adress::registerPython() {
      using namespace espresso::python;

      class_<Adress, shared_ptr<Adress>, bases<Extension> >
        ("integrator_Adress", init<shared_ptr<System> >())
        .def("connect", &Adress::connect)
        .def("disconnect", &Adress::disconnect)
        ;
    }

  }

}
