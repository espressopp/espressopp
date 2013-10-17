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

using namespace std;

namespace espresso {

  namespace integrator {

    using namespace espresso::iterator;

    ExtVirtualParticles::ExtVirtualParticles(shared_ptr<System> system, shared_ptr<CellList> _cl)
        :Extension(system){
        LOG4ESPP_INFO(theLogger, "construct ExtVirtualParticles");
        type = Extension::ExtVirtualParticles;
        vrealCells=_cl;
        onCellListsChanged();
    }


    ExtVirtualParticles::~ExtVirtualParticles() {
      LOG4ESPP_INFO(theLogger, "~ExtVirtualParticles");
      disconnect();
    }

    void ExtVirtualParticles::disconnect(){
        _initForces.disconnect();
        _beforeIntegrate.disconnect();
        _afterIntegrate.disconnect();
    }

    void ExtVirtualParticles::connect() {

    	/*  needed for updating each step */
    	_afterIntegrate = integrator->aftIntP.connect(boost::bind(&ExtVirtualParticles::updateVParticles, this));
    	_afterUpdateGhosts =  integrator->aftUpdGhosts.connect(boost::bind(&ExtVirtualParticles::rebuildVCellLists, this));

    	/* needed for updating on decomposition event */
        _onCellListsChanged = getSystem()->storage->onCellListsChanged.connect(
              boost::bind(&ExtVirtualParticles::onCellListsChanged, this));

        _beforeDecompose = getSystem()->storage->beforeDecompose.connect(0, boost::bind(&ExtVirtualParticles::updateVParticles, this));

        _onParticlesChanged = getSystem()->storage->onParticlesChanged.connect(1, boost::bind(&ExtVirtualParticles::onParticlesChanged, this));


    }

    void ExtVirtualParticles::initRun(){
    }
    void ExtVirtualParticles::initForces(){

    }

    void ExtVirtualParticles::onParticlesChanged(){
    	rebuildVCellLists();
    }


    void ExtVirtualParticles::onCellListsChanged(){
    	System& system = getSystemRef();
		esutil::Error err(system.comm);

		CellList & vcl = *vrealCells;
		// Clean up
		for (CellList::iterator it = vcl.begin(); it != vcl.end(); it++) {
			(*it)->neighborCells.clear();
			(*it)->particles.clear();
			delete (*it);
		}
		for (CellList::iterator it = vghostCells.begin(); it != vghostCells.end(); it++) {
			(*it)->neighborCells.clear();
			(*it)->particles.clear();
			delete (*it);
		}

		vcl.clear();
		vghostCells.clear();
		cellmap.clear();

		const CellList & cl = system.storage->getRealCells();
		for (CellList::const_iterator it = cl.begin(); it != cl.end(); it++) {
			Cell * cellCopy= new Cell();
			vcl.push_back(cellCopy);
			cellmap.insert(std::make_pair<Cell*, Cell*>(*it, cellCopy));
		}
		cout << "Real Cells : " << cl.size() << endl;
		const CellList & gcl = system.storage->getGhostCells();
		for (CellList::const_iterator it = gcl.begin(); it != gcl.end(); it++) {
			Cell * cellCopy= new Cell();
			vghostCells.push_back(cellCopy);
			cellmap.insert(std::make_pair<Cell*, Cell*>(*it, cellCopy));
		}

		std::map<Cell*, Cell*>::iterator cellmap_it;
		for (CellList::const_iterator it = cl.begin(); it != cl.end(); it++) {

			Cell * cellCopy;
			/* Find mapped copy CellList*/
			cellmap_it = cellmap.find(*it);
			if (cellmap_it == cellmap.end()) {
				std::stringstream msg;
				msg << "Missing real cell in ExtVirtualParticle";
				err.setException(msg.str());
				exit(0);
			} else {
				cellCopy = cellmap_it->second;
			}
			NeighborCellList &nb = ((*it)->neighborCells);
			for (NeighborCellList::iterator nbit = nb.begin(); nbit != nb.end();
					nbit++) {
				//lookup nb cells
				NeighborCellInfo & nbinfo = (*nbit);
				//see if we have copied this cell already
				cellmap_it = cellmap.find(nbinfo.cell);
				//NeighborCellInfo * mappednb;
				if (cellmap_it == cellmap.end()) {
					std::stringstream msg;
					msg << "Missing cell in ExtVirtualParticle";
					err.setException(msg.str());
					exit(0);

				} else {
					// store the reference
					cellCopy->neighborCells.push_back(NeighborCellInfo(*cellmap_it->second,nbinfo.useForAllPairs));
				}

			}
		}

    }

	void ExtVirtualParticles::updateVParticles() {
    	System& system = getSystemRef();
		esutil::Error err(system.comm);
		if(cellmap.size()==0){
			cout<< "Cellmap has size 0, rebuild!" << endl;
			onCellListsChanged();
		}

		if (!fixedTupleList){
			std::stringstream msg;
			msg << "FixedTupleList not set in ExtVirtualParticles";
			err.setException( msg.str() );
		}

		std::map<Cell*, Cell*>::iterator cellmap_it;
		const CellList & cl = system.storage->getRealCells();
		for (CellList::const_iterator it = cl.begin(); it != cl.end(); it++) {
			Cell * cellCopy;
			/* Check if this this cell has been copied already in NeighborCellList?*/
			cellmap_it = cellmap.find(*it);
			if (cellmap_it == cellmap.end()) {
				std::stringstream msg;
				msg << "Missing real cell in ExtVirtualParticle";
				err.setException( msg.str() );
			}else{
				cellCopy=cellmap_it->second;
			}
			//clear particles from last built
			//cellCopy->particles.clear();
			ParticleList &part = (*it)->particles;

			for (ParticleList::iterator it2 = part.begin(); it2 != part.end();it2++) {
				Particle &p = (*it2);

				if (std::find(vp_types.begin(), vp_types.end(), p.getType())
						!= vp_types.end()) {
					// This is a virtual particle, update its position based on the COM
					Real3D com =fixedTupleList->calcTupleCOM(p.getId());
					p.setPos(com);

				}
			}
		}

		err.checkException();

	}

void ExtVirtualParticles::rebuildVCellLists() {
	System& system = getSystemRef();
	esutil::Error err(system.comm);

	if (cellmap.size() == 0) {
		cout << "Cellmap has size 0, rebuild!" << endl;
		onCellListsChanged();
	}

	if (!fixedTupleList) {
		std::stringstream msg;
		msg << "FixedTupleList not set in ExtVirtualParticles";
		err.setException(msg.str());
		exit(0);
	}

	std::map<Cell*, Cell*>::iterator cellmap_it;
	const CellList & cl = system.storage->getLocalCells();
	for (CellList::const_iterator it = cl.begin(); it != cl.end(); it++) {
		Cell * cellCopy;
		/* Check if this this cell has been copied already in NeighborCellList?*/
		cellmap_it = cellmap.find(*it);
		if (cellmap_it == cellmap.end()) {
			std::stringstream msg;
			msg << "Missing real cell in ExtVirtualParticle";
			err.setException(msg.str());
		} else {
			cellCopy = cellmap_it->second;
		}
		//clear particles from last built
		cellCopy->particles.clear();
		ParticleList &part = (*it)->particles;

		for (ParticleList::iterator it2 = part.begin(); it2 != part.end();
				it2++) {
			Particle &p = (*it2);

			if (std::find(vp_types.begin(), vp_types.end(), p.getType())
					!= vp_types.end()) {
				cellCopy->particles.push_back(p); //copy of p is made
			}
		}

	}

}



    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void ExtVirtualParticles::registerPython() {
      using namespace espresso::python;

      class_<ExtVirtualParticles, shared_ptr<ExtVirtualParticles>, bases<Extension> >
        ("integrator_ExtVirtualParticles", init<shared_ptr<System>, shared_ptr<CellList>  >())
        .def("connect", &ExtVirtualParticles::connect)
        .def("disconnect", &ExtVirtualParticles::disconnect)
        .def("addVirtualParticleType", &ExtVirtualParticles::addVirtualParticleType)
        .def("setFixedTupleList", &ExtVirtualParticles::setFixedTupleList)
        //.def("getCellList", &ExtVirtualParticles::getCellList, return_value_policy< reference_existing_object >())
        ;
    }

  }

}

