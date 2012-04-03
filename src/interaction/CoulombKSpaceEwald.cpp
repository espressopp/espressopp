#include "python.hpp"
#include <boost/signals2.hpp>
#include "CoulombKSpaceEwald.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class CellListAllParticlesInteractionTemplate <CoulombKSpaceEwald> CellListCoulombKSpaceEwald;

    CoulombKSpaceEwald::CoulombKSpaceEwald(shared_ptr< System > _system, real _prefactor, real _alpha, int _kmax){
      system = _system;
      prefactor = _prefactor;
      alpha = _alpha;
      kmax  = _kmax;
      
      I = Tensor(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
      sum = NULL;
      totsum = NULL;
      
      preset();
      // getParticleNumber(); // geting the number of particles for the current node // it's done in preset
      
      // This function calculates the square of all particle charges. It should be called ones,
      // if the total number of particles doesn't change.
      count_charges(system->storage->getRealCells()); 

      // make a connection to boundary conditions to invoke recalculation of KVec if box dimensions change
      connectionRecalcKVec = system->bc->onBoxDimensionsChanged.connect(boost::bind(&CoulombKSpaceEwald::preset, this));
      // make a connection to storage to get number of particles
      connectionGetParticleNumber = system->storage->onParticlesChanged.connect(boost::bind(&CoulombKSpaceEwald::getParticleNumber, this));
    }
    
    CoulombKSpaceEwald::~CoulombKSpaceEwald(){
      delete [] sum;
      delete [] totsum;
      sum = NULL;
      totsum = NULL;
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void CoulombKSpaceEwald::registerPython() {
      using namespace espresso::python;

      class_< CoulombKSpaceEwald, bases< Potential > >
    	("interaction_CoulombKSpaceEwald", init< shared_ptr< System >, real, real, int >())
    	.add_property("prefactor", &CoulombKSpaceEwald::getPrefactor, &CoulombKSpaceEwald::setPrefactor)
    	.add_property("alpha", &CoulombKSpaceEwald::getAlpha, &CoulombKSpaceEwald::setAlpha)
    	.add_property("kmax", &CoulombKSpaceEwald::getKMax, &CoulombKSpaceEwald::setKMax)
      ;

      class_< CellListCoulombKSpaceEwald, bases< Interaction > >
        ("interaction_CellListCoulombKSpaceEwald",	init< shared_ptr< storage::Storage >, shared_ptr< CoulombKSpaceEwald > >())
        .def("getPotential", &CellListCoulombKSpaceEwald::getPotential)
	  ;

    }
    
  }
}
