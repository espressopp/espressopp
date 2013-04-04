#include "python.hpp"
//#include <boost/signals2.hpp>
#include "CoulombKSpaceP3M.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class CellListAllParticlesInteractionTemplate <CoulombKSpaceP3M> CellListCoulombKSpaceP3M;

    CoulombKSpaceP3M::
    CoulombKSpaceP3M(shared_ptr< System > _system, real _prefactor, real _alpha, int _kmax){
      system = _system;
      
      preset();
      
      // This function calculates the square of all particle charges. It should be called ones,
      // if the total number of particles doesn't change.
      count_charges(system->storage->getRealCells()); 

      /* make a connection to boundary conditions to invoke recalculation of KVec if box
         dimensions change
      */
      connectionRecalcKVec = system->bc->onBoxDimensionsChanged.
              connect(boost::bind(&CoulombKSpaceP3M::preset, this));
      // make a connection to storage in order to get number of particles
      connectionGetParticleNumber = system->storage->onParticlesChanged.
              connect(boost::bind(&CoulombKSpaceP3M::getParticleNumber, this));
    }
    
    CoulombKSpaceP3M::~CoulombKSpaceP3M(){
      /*
      delete [] sum;
      delete [] totsum;
      sum = NULL;
      totsum = NULL;
      */
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void CoulombKSpaceP3M::registerPython() {
      using namespace espresso::python;

      class_< CoulombKSpaceP3M, bases< Potential > >
    	("interaction_CoulombKSpaceP3M", init< shared_ptr< System >, real, real, int >())
    	//.add_property("prefactor", &CoulombKSpaceP3M::getPrefactor, &CoulombKSpaceP3M::setPrefactor)
    	//.add_property("alpha", &CoulombKSpaceP3M::getAlpha, &CoulombKSpaceP3M::setAlpha)
    	//.add_property("kmax", &CoulombKSpaceP3M::getKMax, &CoulombKSpaceP3M::setKMax)
      ;

      class_< CellListCoulombKSpaceP3M, bases< Interaction > >
        ("interaction_CellListCoulombKSpaceP3M",	init< shared_ptr< storage::Storage >, shared_ptr< CoulombKSpaceP3M > >())
        .def("getPotential", &CellListCoulombKSpaceP3M::getPotential)
	  ;

    }
    
  }
}
