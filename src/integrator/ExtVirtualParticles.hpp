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


namespace espresso {

  namespace integrator {

      class ExtVirtualParticles : public Extension {

      public:

        ExtVirtualParticles(shared_ptr<System> system);

        ~ExtVirtualParticles();

        void rebuildVCellLists(); // update local cellLists (fill them with virtual particles and update their center of mass)

        /** Register this class so it can be used from Python. */
        static void registerPython();

        void setFixedTupleList(shared_ptr<FixedTupleList>  ftl){
        	fixedTupleList=ftl;
        }

        void addVirtualParticleType(int type){
        	vp_types.push_back(type);
        }

      private:

        boost::signals2::connection _initForces, _integrate1, _integrate2;
        boost::signals2::connection _onCellListsChanged;

        void integrate1(real&);
        void initForces();
        void integrate2();

        void connect();
        void disconnect();

        void onCellListsChanged(); // called by signal onCellListsChanged, will copy cellList structure to vrealCells vghostCells

        std::vector<int> vp_types; // Types of virtual particles
        CellList vrealCells, vghostCells;
        shared_ptr<FixedTupleList> fixedTupleList;
        std::map<Cell*, Cell*> cellmap; // need to know which one of the new cells corresponds to which original cell
      };

  }

}

#endif	/* EXTVIRTUALPARTICLES_HPP */

