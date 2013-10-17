#include "python.hpp"
#include <sstream>
#include "FixedTupleList.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"



namespace espresso {

  LOG4ESPP_LOGGER(FixedTupleList::theLogger, "FixedTupleList");

  FixedTupleList::FixedTupleList(shared_ptr< storage::Storage > _storage)
    : storage(_storage), globalTuples()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedTupleList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedTupleList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedTupleList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (0, boost::bind(&FixedTupleList::onParticlesChanged, this));

  }

  FixedTupleList::~FixedTupleList() {

    LOG4ESPP_INFO(theLogger, "~FixedTupleList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
    }

    bool FixedTupleList::
    addTuple(boost::python::list& tuple) {

    	longint pidK; // the pid used as key
    	std::vector<longint> pidstmp; //used to extract from tuple;

    	pidK=boost::python::extract<int>(tuple[0]);
    	for (int i = 1; i < boost::python::len(tuple); ++i) {
    		pidstmp.push_back(boost::python::extract<int>(tuple[i]));
    	}
    	globalTuples.insert(make_pair(pidK, pidstmp));

        /*bool returnVal = true;
        System& system = storage->getSystemRef();
        esutil::Error err(system.comm);
        
        Particle* vp, *at;
        longint pidK; // the pid used as key
        std::vector<Particle*> tmp; // temporary vector
        std::vector<longint> pids; //used to extract from tuple;
        std::vector<longint> pidstmp; // temporary vector
        

        for (int i = 0; i < boost::python::len(tuple); ++i) {
            pids.push_back(boost::python::extract<int>(tuple[i]));
        }
        
        tuple::iterator it = pids.begin();
        vp = storage->lookupRealParticle(*it);
        if (!vp) {// Particle does not exist here, return false
            returnVal = false;
        } else {
            pidK = *it; // first pid is key
            for (++it; it != pids.end(); ++it) {
                at = storage->lookupLocalParticle(*it);
                if (!at) { // Particle does not exist here, return false
                    std::stringstream msg;
                    msg << "ERROR: particle " << *it << " not found \n";
                    err.setException(msg.str());
                    returnVal = false;
                    break;
                }
                tmp.push_back(at);
                pidstmp.push_back(*it); // pidK is not in this vector
            }
        }
        err.checkException();

        if (returnVal) {
            this->add(vp, tmp); // add to TupleList
            globalTuples.insert(make_pair(pidK, pidstmp));
        }
        LOG4ESPP_INFO(theLogger, "added fixed tuple to global tuples");

        tmp.clear();
        pids.clear();
        pidstmp.clear();


        return returnVal;*/
    }

    python::list FixedTupleList::getTuples() {
        python::list tuple;
        python::list alltuples;
        for (GlobalTuples::iterator it = globalTuples.begin();
                it != globalTuples.end(); it++) {
            python::list tuple;
            tuple.append((*it).first); // key is also part of the tuple!
            for (tuple::iterator it2=(*it).second.begin(); it2 !=(*it).second.end(); it2++){
                tuple.append(*it2);
            }
            alltuples.append(tuple);
        }

        return alltuples;
    }

   std::vector<Particle *> FixedTupleList::getTupleByID(int id){
	   // TODO: enable fast access through an unordererd map id-> tuple
	   Particle * vp = storage->lookupLocalParticle(id);
	   TupleList::iterator mapit= this->find(vp);
	   if (mapit == this->end()){
		   std::cout << "ERROR: VP with id " << id << "not found." << "local particle id" << vp->id() << std::endl;
		   exit(0);
	   }

	   return mapit->second;
   }

  void FixedTupleList::
    beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
        /*std::vector<longint> toSend;
        longint vpcount=0; // number of virtual particles
        // loop over the particle list
        for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
            longint pidK = pit->id();
            LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pidK << ", find tuples");

            // find particle that has this particle id
            GlobalTuples::const_iterator it = globalTuples.find(pidK);
            if (it != globalTuples.end()) {
            	vpcount++;
                // write the size of the vector
                longint tuplesize = it->second.size();

                int extend=2+tuplesize;
                toSend.reserve(toSend.size()+extend);

                toSend.push_back(pidK);
                std::cout << "send VP" << pidK << std::endl;
                toSend.push_back(tuplesize);
                // how many tuple ids are following
                // iterate through vector and add pids
                for (tuple::const_iterator it2 = it->second.begin();
                        it2 != it->second.end(); ++it2) {
                    toSend.push_back(*it2);
                }
                // delete this pid from the global list
                globalTuples.erase(pidK);

            }
        }

        buf.write(vpcount);

        for(std::vector<longint>::iterator it=toSend.begin();it!=toSend.end();it++){
        	//manual writing to avoid writing the size of the vector (not needed)
        	buf.write(*it);
        }*/

    }

    void FixedTupleList::
    afterRecvParticles(ParticleList &pl, InBuffer& buf) {
        /*std::vector<longint> pids;
        longint tuplesize=0;
        longint vpcount=0;
        longint pidK, id;

        GlobalTuples::iterator it = globalTuples.begin();

        buf.read(vpcount);
        for (int i=0; i<vpcount;i++) {
            buf.read(pidK);
            buf.read(tuplesize);
            pids.reserve(tuplesize);
            for (int j=0; j<tuplesize; j++){
            	buf.read(id);
            	pids.push_back(id);
            }

            it = globalTuples.insert(it, std::make_pair(pidK, pids));
            pids.clear();
        }*/

    }
   Real3D FixedTupleList::calcTupleCOM(int id){
	   System& system = storage->getSystemRef();
	   GlobalTuples::const_iterator it;
	   // TODO: instead of global tuples we could search in  some local list. Local list has to be generated from reals+ghosts
	   it=globalTuples.find(id);
	   Real3D com=Real3D(0,0,0);
	   const bc::BC& bc = *system.bc;  // boundary conditions
	   real m=0.0;

	   if (it!=globalTuples.end()){
	  		   Particle* p;
	  		   real M=0.0;
	  		   tuple::const_iterator it2 = it->second.begin();
	  		   p = storage->lookupLocalParticle(*it2);

	  		   if (id==93) std::cout << "VP " << p->position() << std::endl;
	  		   if (p){
	  			 com=p->getPos();
	  			 system.bc->foldPosition(com);
	  			 M=p->getMass();
	  			it2++;
	  		   }
	  		   for (;it2 != it->second.end(); ++it2) {

	  			 p = storage->lookupLocalParticle(*it2);
	  			if (id==93){
	  				Particle * p2 = storage->lookupRealParticle(*it2);
	  				if (p2) std::cout<< "PPP rank:"<< system.comm->rank() << " I also have to offer p2 id " << *it2  << " "<< p2->getPos() << std::endl;
	  			}
	  			 if (id==93) std::cout << "pos " << p->position() << " comi "<<  com << std::endl;
	  			 Real3D dx;
	  			 bc.getMinimumImageVectorBox(dx, p->getPos(), com);
	  			 m=p->getMass();
	  			 com=com+(m/(M+m))*dx;
	  			 M+=m;
	  		   }
	  		 system.bc->foldPosition(com);
	  		if (id==93) std::cout << "final folded com " << com << std::endl;
	   }else{
		   std::stringstream msg;
		   msg << " tuple with key id " << id << " does not exists on node " << mpiWorld->rank();
		   return com;
	   }



	   /*if (it!=globalTuples.end()){
		   Particle* p;
		   real M=0.0;
		   for (tuple::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
			   p = storage->lookupLocalParticle(*it2);
			   Real3D upos=p->getPos();
			   Int3D img=p->getImageBox();
			   if (id==2919){
				   	   	   	   std::cout << "("<< p->getId() <<  ") ghost:" << p->ghost() << std::endl;
			   				   std::cout << "("<< p->getId() <<  ")" <<" pos " << upos << " img " << img << std::endl;
			   			   }
			   system.bc->unfoldPosition(upos, img);
			   M+=p->getMass();
			   com+=upos*p->getMass();
			   if (id==2919){
				   std::cout << "("<< p->getId() <<  ")" << "upos " << upos << " img " << img << std::endl;
			   }
		   }

		   if (M>0) com/=M;
		   if (id==2919){
		   std::cout << "com" << com << std::endl;
		   }
		   system.bc->foldPosition(com);
		   if (id==2919){
		   		   std::cout << "folded com" << com << std::endl;
		   		   }
		   return com;
	   }
	   else{
		   std::stringstream msg;
		   msg << " tuple with key id " << id << " does not exists on node " << mpiWorld->rank();
		   return com;
	   }*/
    }

   void FixedTupleList::unwrapMinimumImage(int id){
	   System& system = storage->getSystemRef();
	   GlobalTuples::const_iterator it;

	   const bc::BC& bc = *system.bc;  // boundary conditions

	   it=globalTuples.find(id);
	   tuple::const_iterator it2 = it->second.begin();
	   for (;it2 != it->second.end(); ++it2) {

	   }

   }


  void FixedTupleList::onParticlesChanged() {
         // TODO errors should be thrown in a nice way
        LOG4ESPP_INFO(theLogger, "rebuild local particle list from global tuples\n");
        System& system = storage->getSystemRef();
        esutil::Error err(system.comm);
        
        this->clear();
        Particle* vp, * at;
        std::vector<Particle*> tmp;

        GlobalTuples::const_iterator it = globalTuples.begin();

        for (;it != globalTuples.end(); ++it) {
            vp = storage->lookupLocalParticle(it->first);
            if (vp == NULL) continue;
            /*if (vp == NULL) {
            	std::stringstream msg;
                msg << " particle in tuple " << it->first << " does not exists here";
                err.setException( msg.str() );
            }*/


            // iterate through vector in map
            for (tuple::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                at = storage->lookupLocalParticle(*it2);
                if (at == NULL) {
                	std::stringstream msg;
                        msg << " particle in tuple " << *it2 << " does not exists here";
                        err.setException( msg.str() );
                }
                tmp.push_back(at);
            }
            this->add(vp, tmp);
            tmp.clear();
        }
        LOG4ESPP_INFO(theLogger, "regenerated local fixed list from global tuples");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTupleList::registerPython() {

    using namespace espresso::python;

    class_< FixedTupleList, shared_ptr< FixedTupleList > >
      ("FixedTupleList", init< shared_ptr< storage::Storage > >())
      .def("addTuple", &FixedTupleList::addTuple)
     .def("getTuples", &FixedTupleList::getTuples)
      .def("size", &FixedTupleList::size)
     ;
  }
}
