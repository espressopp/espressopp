#include "python.hpp"
#include "ExtVirtualParticles.hpp"
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

    E::ExtVirtualParticles(shared_ptr<System> system)
        :Extension(system){
        LOG4ESPP_INFO(theLogger, "construct ExtVirtualParticles");
        type = Extension::ExtVirtualParticles;
    }


    ExtVirtualParticles::~ExtVirtualParticles() {
      LOG4ESPP_INFO(theLogger, "~ExtVirtualParticles");
      disconnect();
    }

    void ExtVirtualParticles::disconnect(){
        _initForces.disconnect();
        _integrate1.disconnect();
        _integrate2.disconnect();
    }

    void ExtVirtualParticles::connect() {

        // connection to after initForces()
        _initForces = integrator->aftInitF.connect(
                boost::bind(&ExtVirtualParticles::initForces, this));

        // connection to inside of integrate1()
        _integrate1 = integrator->inIntP.connect(
                boost::bind(&ExtVirtualParticles::integrate1, this, _1));

        // connection to after integrate2()
        _integrate2 = integrator->aftIntV.connect(
                boost::bind(&ExtVirtualParticles::integrate2, this));
    }


    void ExtVirtualParticles::initForces(){

    }



    void ExtVirtualParticles::integrate1(real& maxSqDist){

    }


    void ExtVirtualParticles::integrate2() {
        CellList cl= system.storage->getLocalCells();
        CellList::iterator it;
        
        CellList vcl;
        std::map<Cell*, Cell*> cellmap; // need to know which one of the new cells corresponds to which cell
        
        for (;it!=cl.end();cl++){ 
            shared_ptr<Cell> newCell;
            
            NeighborCellList * nb=it->neighborCells;
            for (NeighborCellList::iterator nbit=){
                //create or lookup cell
                
            }
        }
        
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void ExtVirtualParticles::registerPython() {
      using namespace espresso::python;

      class_<ExtVirtualParticles, shared_ptr<ExtVirtualParticles>, bases<Extension> >
        ("integrator_ExtVirtualParticles", init<shared_ptr<System> >())
        .def("connect", &ExtVirtualParticles::connect)
        .def("disconnect", &ExtVirtualParticles::disconnect)
        ;
    }

  }

}

