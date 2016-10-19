/* 
 * File:   ExtVirtualParticles.hpp
 * Author: fritsch
 *
 * Created on August 21, 2013, 11:43 AM
 */

#ifndef EXTVIRTUALPARTICLES_HPP
#define	EXTVIRTUALPARTICLES_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "Cell.hpp"
#include "boost/signals2.hpp"
#include "FixedTupleList.hpp"


namespace espressopp {

  namespace integrator {



      class ExtVirtualParticles : public Extension {
    	  /* Provides an extension class which takes care of updating 'virtual' particles in each step.
    	   * Currently, the virtual particle is placed in the center of mass of the constituting particles.
    	   * A cell list is provided which contains the virtual particles. This can be used to build verlet lists using the VirtualVerletList
    	   * class.
    	   */

      public:

        ExtVirtualParticles(shared_ptr<System> system, shared_ptr<CellList> _cl);

        ~ExtVirtualParticles();



        /** Register this class so it can be used from Python. */
        static void registerPython();

        void setFixedTupleList(shared_ptr<FixedTupleList>  ftl){
        	fixedTupleList=ftl;
        }

        void addVirtualParticleType(int type){
        	vp_types.push_back(type);
        }

      private:

        boost::signals2::connection _initForces, _runInit, _beforeIntegrate, _afterIntegrate, _afterUpdateGhosts, _beforeDecompose;
        boost::signals2::connection _onCellListsChanged, _onParticlesChanged;

        void initForces();

        void connect();
        void disconnect();

        void initRun();

        void onCellListsChanged(); // called by signal onCellListsChanged, will copy cellList structure to vrealCells vghostCells
        void rebuildVCellLists(); // update local cellLists (fill them with virtual particles)
        void updateVParticles(); /*update the positions of the virtual particles */
        void onParticlesChanged();

        std::vector<int> vp_types; // Types of virtual particles
        CellList vghostCells;
        shared_ptr<CellList> vrealCells;
        shared_ptr<FixedTupleList> fixedTupleList;
        std::map<Cell*, Cell*> cellmap; // need to know which one of the new cells corresponds to which original cell
      };

  }

}

#endif	/* EXTVIRTUALPARTICLES_HPP */

