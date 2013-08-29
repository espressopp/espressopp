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

        _onCellListsChanged = getSystem()->storage->onCellListsChanged.connect(
              boost::bind(&ExtVirtualParticles::onCellListsChanged, this));
    }


    void ExtVirtualParticles::initForces(){

    }



    void ExtVirtualParticles::integrate1(real& maxSqDist){

    }


    void ExtVirtualParticles::integrate2() {
        rebuildVCellLists();

    }

    void ExtVirtualParticles::onCellListsChanged(){
    	System& system = getSystemRef();
		esutil::Error err(system.comm);
		cout << "CellListsChanged" << endl;

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
			} else {
				cellCopy = cellmap_it->second;
			}
			NeighborCellList &nb = ((*it)->neighborCells);
			for (NeighborCellList::iterator nbit = nb.begin(); nbit != nb.end();
					nbit++) {
				//create or lookup cells
				NeighborCellInfo & nbinfo = (*nbit);
				//see if we have copied this cell already
				cellmap_it = cellmap.find(nbinfo.cell);
				//NeighborCellInfo * mappednb;
				if (cellmap_it == cellmap.end()) {
					std::stringstream msg;
					msg << "Missing cell in ExtVirtualParticle";
					err.setException(msg.str());
					// not copied yet, make an emtpy copy and store it
					//Cell* mapped_cell = new Cell();
					//mappednb= new NeighborCellInfo(mapped_cell,nbinfo.useForAllPairs);
					//cellmap.insert(std::make_pair<Cell*, Cell*>(nbinfo.cell, mapped_cell));

				} else {
					// store the reference
					//mappednb = new NeighborCellInfo(*cellmap_it->second, nbinfo.useForAllPairs);
					cellCopy->neighborCells.push_back(
							NeighborCellInfo(*cellmap_it->second,
									nbinfo.useForAllPairs));
				}

			}
		}

		cout << "lenth cellmap" << cellmap.size() << endl;
    }

	void ExtVirtualParticles::rebuildVCellLists() {
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
				err.setException( msg.str() );
			}else{
				cellCopy=cellmap_it->second;
			}
			//clear particles from last built
			cellCopy->particles.clear();
			ParticleList &part = (*it)->particles;

			for (ParticleList::iterator it2 = part.begin(); it2 != part.end();it2++) {
				std::vector<int>::iterator idit;
				Particle &p = (*it2);

				if (std::find(vp_types.begin(), vp_types.end(), p.getType())
						!= vp_types.end()) {
					// This is a virtual particle, update its position based on the COM
					Real3D com =fixedTupleList->calcTupleCOM(p.getId());
					p.setPos(com);
					cellCopy->particles.push_back(p);
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

