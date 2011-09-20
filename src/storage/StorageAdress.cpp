#include "python.hpp"

//#include <algorithm>

#include "log4espp.hpp"

#include "System.hpp"
#include "StorageAdress.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListIterator.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "Particle.hpp"
#include "Buffer.hpp"

using namespace boost;
using namespace espresso::iterator;

namespace espresso {
 
  namespace storage {
    LOG4ESPP_LOGGER(StorageAdress::logger, "StorageAdress");


    StorageAdress::StorageAdress(shared_ptr< System> system)
    : Storage(system)
    {
      //logger.setLevel(log4espp::Logger::TRACE);
      LOG4ESPP_INFO(logger, "Created new storage object for a system, has buffers");
    }

    StorageAdress::~StorageAdress() {}


    void StorageAdress::invalidateGhosts() {
      for(CellListIterator it(getGhostCells());
	  it.isValid(); ++it) {
	/* remove only ghosts from the hash if the localParticles hash
	   actually points to the ghost.  If there are local ghost cells
	   to implement pbc, the real particle will be the one accessible
	   via localParticles.
	*/
          removeFromLocalParticles(&(*it), true);
      }

      // for AdResS
      //std::cout << getSystem()->comm->rank() << ": ----- CLEAR GHOST ADRESS PARTICLES (AdrATParticlesG) ----\n";
      //AdrATParticlesG.clear();
      clearAdrATParticlesG();
      // clear ghost tuples
      FixedTupleList::iterator it = fixedtupleList->begin();
      for (;it != fixedtupleList->end(); ++it) {
    	  Particle* vp = it->first;
    	  if (vp->ghost()) {
    		  //std::cout << "erasing ghost particle in tuple: " << vp->id() << "-" << vp->ghost() << "\n";
    		  fixedtupleList->erase(it);
    	  }
      }
    }

    void StorageAdress::decompose() {
      //std::cout << " ---- decompose ----\n";
      invalidateGhosts();
      decomposeRealParticles();
      //std::cout << getSystem()->comm->rank() << ": (onTuplesChanged) ";
      onTuplesChanged(); // for AdResS, renamed to not confuse with bonds
      //std::cout << " ---- exchange ghosts ---- \n";
      exchangeGhosts();
      //std::cout << getSystem()->comm->rank() << ": ";
      onParticlesChanged();
    }

    void StorageAdress::packPositionsEtc(OutBuffer &buf,
				   Cell &_reals, int extradata, const Real3D& shift) {
      ParticleList &reals  = _reals.particles;

      for(ParticleList::iterator src = reals.begin(), end = reals.end(); src != end; ++src) {

        buf.write(*src, extradata, shift);


        /*std::cout << getSystem()->comm->rank() << ": -> " << src->id() << "-" << src->ghost() <<
        		" (" << src->getPos() << ") type " << src->getType() << " @ " << &(*src) << " \n";*/


        // for AdResS
		// write all AT particles belonging to this VP into the buffer
        FixedTupleList::iterator it;
		it = fixedtupleList->find(&(*src));
		if (it != fixedtupleList->end()) {
			std::vector<Particle*> atList;
			atList = it->second;

			int size = atList.size();
			buf.write(size); // write size of vector first

			for (std::vector<Particle*>::iterator itv = atList.begin();
				  itv != atList.end(); ++itv) {
				Particle &at = **itv;

				/*std::cout << getSystem()->comm->rank() << ":  ---> " << at.id() << "-" << at.ghost() <<
						" (" << at.position() << ") type " << at.type() << " \n";*/

				/*
				if (at.type() == 0) { // just for testing
					std::cout << "SERIOUS ERROR, particle is of wrong type\n";
					exit(1);
					return;
				}*/

				buf.write(at, 1, shift); // we force extradata to 1
			}
		}
		else {
			std::cout << getSystem()->comm->rank() << ": packposetc " << "VP particle "<< src->id() << "-" << src->ghost() << " not found in tuples!\n";
			exit(1);
			return;
		}

      }
    }

    void StorageAdress::unpackPositionsEtc(Cell &_ghosts, InBuffer &buf, int extradata) {
      ParticleList &ghosts  = _ghosts.particles;

      for(ParticleList::iterator dst = ghosts.begin(), end = ghosts.end(); dst != end; ++dst) {

    	//std::cout << getSystem()->comm->rank() << ": buf.read(particle, extradata) (unpackPosEtc) \n";
        buf.read(*dst, extradata);

        if (extradata & DATA_PROPERTIES) {
        	updateInLocalParticles(&(*dst), true);
        }

        dst->ghost() = 1;


        // for AdResS

        /*std::cout << getSystem()->comm->rank() << ": <- " << dst->id() << "-" << dst->ghost() <<
        		" (" << dst->getPos() << ") type " << dst->getType() << "\n";*/

        //int numAT = fixedtupleList->getNumPart(dst->id()); // does not work, not sure why
        int numAT;
        buf.read(numAT); // read number of AT particles

        FixedTupleList::iterator it;
        it = fixedtupleList->find(&(*dst));
        if (it != fixedtupleList->end()) {
            std::vector<Particle*> atList = it->second;

            for (std::vector<Particle*>::iterator itv = atList.begin();
                    itv != atList.end(); ++itv) {
                Particle &atg = **itv;
                buf.read(atg, 1); // we force extradata to 1 (although not necessary here...)

                /*std::cout << getSystem()->comm->rank() << ":  <--- " << atg.id() << "-" << atg.ghost() <<
                                " (" << atg.position() << ") type " << atg.type() <<
                                " extradata: " << extradata << "\n";*/
            }
        }
        else {
            std::vector<Particle*> tmp;
            Particle* atg; // atom, ghost
            Particle tmpatg; // temporary particle, to be inserted into adr. at. ghost part.

            for (int i = 1; i <= numAT; ++i) {
                // read the AT partcles
                //std::cout << getSystem()->comm->rank() << ": buf.read(AT particle, extradata) (unpackPosEtc), i: " << i << "\n";
                buf.read(tmpatg, 1); // we force extradata to 1

                atg = appendUnindexedAdrParticle(getAdrATParticlesG(), tmpatg);
                atg->ghost() = 1;

                tmp.push_back(atg);

                /*std::cout << getSystem()->comm->rank() << ":  <--- " << atg->id() << "-" << atg->ghost() <<
                                    " (" << atg->getPos() << ") type " << atg->getType() <<
                                    " extradata: " << extradata << "\n";*/

                /*if (atg->getType() == dst->getType()) { // for testing purposes
                    std::cout << "SERIOUS ERROR, particle is of wrong type\n";
                    exit(-1);
                    return;
                }*/
            }
            fixedtupleList->insert(std::make_pair(&(*dst), tmp));
            tmp.clear();
        }




      }
    }

    void StorageAdress::copyRealsToGhosts(Cell &_reals, Cell &_ghosts,
				    int extradata,
				    const Real3D& shift)
    {
      ParticleList &reals  = _reals.particles;
      ParticleList &ghosts = _ghosts.particles;

      ghosts.resize(reals.size());

      //std::cout << "Copy reals to ghosts ... \n";

      for(ParticleList::iterator src = reals.begin(), end = reals.end(), dst = ghosts.begin();
              src != end; ++src, ++dst) {
        dst->copyAsGhost(*src, extradata, shift);

        // for AdResS
        copyGhostTuples(*src, *dst, extradata, shift);
      }
    }


    void StorageAdress::copyGhostTuples(Particle& src, Particle& dst, int extradata, const Real3D& shift) {

        // create ghosts of particles in tuples
        FixedTupleList::iterator it;
        it = fixedtupleList->find(&src);
        if (it != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it->second;

            Particle* atg; // atom, ghost
            std::vector<Particle*> tmp; // temporary vector

            for (std::vector<Particle*>::iterator itv = atList.begin();
                  itv != atList.end(); ++itv) {
                Particle &at = **itv;

                Particle n; // temporary particle, to be inserted into adr. at. ghost part.
                n.id() = at.getId();
                n.type() = at.getType();


                // see whether the array was resized; STL hack
                //Particle *begin = &AdrATParticlesG.front();

                atg = appendUnindexedAdrParticle(getAdrATParticlesG(), n);
                atg->copyAsGhost(at, extradata, shift);

                tmp.push_back(atg);

                /*if (begin != &AdrATParticlesG.front())
                    std::cout << "\n  -- AdrATParticlesG array resized!! --\n\n";*/

            }
            fixedtupleList->insert(std::make_pair(&dst, tmp));
            tmp.clear();
        }
        else {
            std::cout << "copyGhostTuples: VP particle "<< src.id() << "-" << src.ghost() << " (" << src.position() << ")" << " not found in tuples!\n";
            exit(1);
            return;
        }
    }

    /* -- this is now solved in FixedTupleList.cpp in onparticleschanged()
    void StorageAdress::foldAdrPartCoor(Particle& part, Real3D& oldpos, int coord) {

        if (part.position()[coord] != oldpos[coord]) {
            real moved = oldpos[coord] - part.position()[coord];

            FixedTupleList::iterator it;
            it = fixedtupleList->find(&part);
            if (it != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it->second;

                for (std::vector<Particle*>::iterator itv = atList.begin();
                      itv != atList.end(); ++itv) {
                    Particle &at = **itv;

                    //std::cout << " updating position for AT part " << at.id() << " (" << at.position() << ") ";

                    at.position()[coord] = at.position()[coord] - moved;
                    //getSystem()->bc->foldCoordinate(at.position(), at.image(), coord);

                    //std::cout << "to (" << at.position() << ")\n";
                }
            }
            else {
                std::cout << getSystem()->comm->rank() << ": foldAdrPartCoor "
                        << "VP particle "<< part.id() << "-" << part.ghost() << " not found in tuples!\n";
                exit(1);
                return;
            }
        }
    }
    */

    void StorageAdress::packForces(OutBuffer &buf, Cell &_ghosts) {

      ParticleList &ghosts = _ghosts.particles;
  
      for(ParticleList::iterator src = ghosts.begin(), end = ghosts.end(); src != end; ++src) {

        buf.write(src->particleForce());

        LOG4ESPP_TRACE(logger, "from particle " << src->id() << ": packing force " << src->force());


        /*std::cout << getSystem()->comm->rank() << ": " << " from particle " << src->id()
        		<< " packing force " << src->force() << "\n";*/



        // for AdResS
		// write all AT particle forces belonging to this VP into the buffer
        FixedTupleList::iterator it;
		it = fixedtupleList->find(&(*src));
		if (it != fixedtupleList->end()) {
			std::vector<Particle*> atList;
			atList = it->second;

			for (std::vector<Particle*>::iterator itv = atList.begin();
				  itv != atList.end(); ++itv) {
				Particle &at = **itv;

				/*std::cout << getSystem()->comm->rank() << ":  ---> " << at.id() << "-" << at.ghost() <<
						" (" << at.position() << ") type " << at.type() << " \n";*/

				buf.write(at.particleForce());
			}
		}
		else {
			std::cout << getSystem()->comm->rank() << ": packforces " << "VP particle "<< src->id() << "-" << src->ghost() << " not found in tuples!\n";
			exit(1);
			return;
		}
      }
    }

    void StorageAdress::unpackAndAddForces(Cell &_reals, InBuffer &buf)
    {
      LOG4ESPP_DEBUG(logger, "add forces from buffer to cell "
		     << (&_reals - getFirstCell()));

      ParticleList &reals = _reals.particles;

      for(ParticleList::iterator dst = reals.begin(), end = reals.end(); dst != end; ++dst) {
    	  ParticleForce f;
    	  //std::cout << getSystem()->comm->rank() << ": buf.read(force) (unpackAndAddForces) \n";
    	  buf.read(f);
    	  LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": unpacking force "
		       << f.f() << " and adding to " << dst->force());
    	  dst->particleForce() += f;



    	  /*std::cout << getSystem()->comm->rank() << ": " << " for particle " << dst->id() << " unpacking force "
    		  << " and adding to " << dst->force() << "\n";*/


		  // for AdResS

    	  // iterate through atomistic particles in fixedtuplelist
		  FixedTupleList::iterator it;
		  it = fixedtupleList->find(&(*dst));

		  if (it != fixedtupleList->end()) {

			 std::vector<Particle*> atList1;
			 atList1 = it->second;

			 //std::cout << "AT forces ...\n";
			 for (std::vector<Particle*>::iterator itv = atList1.begin(); itv != atList1.end(); ++itv) {
				 Particle &p3 = **itv;
				 //std::cout << getSystem()->comm->rank() << ": buf.read(AT force) (unpackAndAddForces) \n";
				 buf.read(f);
				 p3.particleForce() += f;
			 }
		  }
		  else {
			 std::cout << " unpackForces: one of the VP particles not found in tuples: " << dst->id() << "-" << dst->ghost();
			 exit(1);
			 return;
		  }

      }
    }

    void StorageAdress::addGhostForcesToReals(Cell &_ghosts, Cell &_reals)
    {

      ParticleList &reals  = _reals.particles;
      ParticleList &ghosts = _ghosts.particles;

      for(ParticleList::iterator dst = reals.begin(), end = reals.end(), src = ghosts.begin();
              dst != end; ++dst, ++src) {
          LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": adding force "
		       << src->force() << " to " << dst->force());

          dst->particleForce() += src->particleForce();

          // for AdResS
          addAdrGhostForcesToReals(*src, *dst);
      }
    }



    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    StorageAdress::registerPython() {
      using namespace espresso::python;


      class_< StorageAdress, boost::noncopyable >("storage_StorageAdress", no_init)

	.def("addParticle", &StorageAdress::addParticle,
	     return_value_policy< reference_existing_object >())

    .def("addAdrATParticle", &StorageAdress::addAdrATParticle,
         return_value_policy< reference_existing_object >())

    .def("setFixedTuples", &StorageAdress::setFixedTuples)

	.def("lookupLocalParticle", &StorageAdress::lookupLocalParticle,
	     return_value_policy< reference_existing_object >())

	.def("lookupRealParticle", &StorageAdress::lookupRealParticle,
	     return_value_policy< reference_existing_object >())

	.def("decompose", &StorageAdress::decompose)
        .add_property("system", &StorageAdress::getSystem)
	;
    }

  }
}
