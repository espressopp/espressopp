#include "python.hpp"
#include <boost/signals2.hpp>
#include "EwaldKSpace.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class CellListAllParticlesInteractionTemplate <EwaldKSpace> CellListEwaldKSpace;

    EwaldKSpace::EwaldKSpace(shared_ptr< System > _system, real _prefactor, real _alpha, int _kmax){
      system = _system;
      prefactor = _prefactor;
      alpha = _alpha;
      kmax  = _kmax;
      preset();
      // getParticleNumber(); // geting the number of particles for the current node // it's done in preset
      
      // This function calculates the square of all particle charges. It should be called ones,
      // if the total number of particles doesn't change.
      count_charges(system->storage->getRealCells()); 

      // make a connection to boundary conditions to invoke recalculation of KVec if box dimensions change
      connectionRecalcKVec = system->bc->onBoxDimensionsChanged.connect(boost::bind(&EwaldKSpace::preset, this));
      // make a connection to storage to get number of particles
      connectionGetParticleNumber = system->storage->onParticlesChanged.connect(boost::bind(&EwaldKSpace::getParticleNumber, this));
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void EwaldKSpace::registerPython() {
      using namespace espresso::python;

      class_< EwaldKSpace, bases< Potential > >
    	("interaction_EwaldKSpace", init< shared_ptr< System >, real, real, int >())
    	.add_property("prefactor", &EwaldKSpace::getPrefactor, &EwaldKSpace::setPrefactor)
    	.add_property("alpha", &EwaldKSpace::getAlpha, &EwaldKSpace::setAlpha)
    	.add_property("kmax", &EwaldKSpace::getKMax, &EwaldKSpace::setKMax)
      ;

      class_< CellListEwaldKSpace, bases< Interaction > >
        ("interaction_CellListEwaldKSpace",	init< shared_ptr< storage::Storage >, shared_ptr< EwaldKSpace > >())
        .def("getPotential", &CellListEwaldKSpace::getPotential)
	  ;

    }
    
  }
}
