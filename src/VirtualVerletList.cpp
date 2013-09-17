#include "python.hpp"
#include "VirtualVerletList.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "esutil/Error.hpp"


namespace espresso {

  using namespace espresso::iterator;

  LOG4ESPP_LOGGER(VirtualVerletList::theLogger, "VirtualVerletList");

/*-------------------------------------------------------------*/

  // cut is a cutoff (without skin)
  VirtualVerletList::VirtualVerletList(shared_ptr<System> system, real _cut,
		  shared_ptr<FixedTupleList> _ftpl, bool rebuildVL) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VirtualVerletList, cut = " << _cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;
    ftpl = _ftpl;

    //if (rebuildVL && cellList ) rebuild(); // not called if exclutions are provided

  
    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VirtualVerletList::rebuild, this));

    vlArray = esutil::Array2D<shared_ptr<VerletList>, esutil::enlarge>(0, 0, shared_ptr<VerletList>());
  }
  
  real VirtualVerletList::getVerletCutoff(){
    return cutVerlet;
  }
  
  void VirtualVerletList::connect()
  {

  // make a connection to System to invoke rebuild on resort
  connectionResort = getSystem()->storage->onParticlesChanged.connect(
      boost::bind(&VirtualVerletList::rebuild, this));
  }

  void VirtualVerletList::disconnect()
  {

  // disconnect from System to avoid rebuild on resort
  connectionResort.disconnect();
  }

  /*-------------------------------------------------------------*/
  
  void VirtualVerletList::rebuild()
  {
  	System& system = getSystemRef();
	esutil::Error err(system.comm);

	/*if(!cellList){
		std::stringstream msg;
		msg << "VirtualVerleList: no cell list set!";
		err.setException( msg.str() );
		exit(0);
	}*/

    //real cutVerlet = cut + getSystem() -> getSkin();
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    
    vlPairs.clear();

    //CellList cl = getSystem()->storage->getRealCells();
    CellList &cl = (*cellList);
    std::cout << "VVL rebuild: " << cl.size() << std::endl;

    LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      checkPair(*it->first, *it->second);

      LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " and " << it->second->id());
    }
    
    builds++;
    LOG4ESPP_DEBUG(theLogger, "rebuilt VirtualVerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlPairs.size());
  }
  

  /*-------------------------------------------------------------*/
  
  void VirtualVerletList::checkPair(Particle& pt1, Particle& pt2)
  {

    Real3D d = pt1.position() - pt2.position();
    real distsq = d.sqr();

    //std::cout << "Handling vp pair" << pt1.id() << " " << pt2.id() << std::endl;

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id()
                   << " @ " << pt1.position() 
		   << " - p2: " << pt2.id() << " @ " << pt2.position()
		   << " -> distsq = " << distsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (both directions)
    if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;


    //std::cout << "sizex " << vlArray.size_x() << "sizey " << vlArray.size_y() << std::endl;
    //std::cout << "t1 " << pt1.type() << "t2 " << pt2.type() << std::endl;
    FixedTupleList::iterator it;

    std::cout << "Unpacking pair ids " << pt1.id() << "  " << pt2.id() << std::endl;
    std::vector<Particle *> tup1=ftpl->getTupleByID(pt1.id());
    std::vector<Particle *> tup2=ftpl->getTupleByID(pt2.id());


    std::vector<Particle*>::iterator pit1= tup1.begin();


    for(; pit1!=tup1.end(); pit1++){
    	for(std::vector<Particle*>::iterator pit2= tup2.begin(); pit2!=tup2.end(); pit2++){
    		Particle* p1=*pit1;
    		Particle* p2=*pit2;
    		std::cout<< "Mapping t1" << p1->type() << " t2 " <<  p2->type() << " id1 " << p1->id() << " id2 " << p2->id();

    		if (p1->type()< vlArray.size_x() && p2->type()< vlArray.size_y()){
    				shared_ptr<VerletList> mappedVl=vlArray.at(p1->type(), p2->type());
    				if (mappedVl){
    					std::cout<<" to " << mappedVl.get();
    					mappedVl->vlPairs.add(p1, p2);
    				}
    				std::cout << std::endl;
    		}
    	}
    }


    //vlPairs.add(pt1, pt2); // add pair to Verlet List
  }
  
  /*-------------------------------------------------------------*/
  
  int VirtualVerletList::totalSize() const
  {
    System& system = getSystemRef();
    int size = localSize();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
  }

  int VirtualVerletList::localSize() const
  {
    System& system = getSystemRef();
    return vlPairs.size();
  }

  python::tuple VirtualVerletList::getPair(int i) {
	  if (i <= 0 || i > vlPairs.size()) {
	    std::cout << "ERROR VirtualVerletList pair " << i << " does not exists" << std::endl;
	    return python::make_tuple();
	  } else {
	    return python::make_tuple(vlPairs[i-1].first->id(), vlPairs[i-1].second->id());
	  }
  }


  bool VirtualVerletList::exclude(longint pid1, longint pid2) {

      exList.insert(std::make_pair(pid1, pid2));

      return true;
  }
  

  /*-------------------------------------------------------------*/
  
  VirtualVerletList::~VirtualVerletList()
  {
    LOG4ESPP_INFO(theLogger, "~VirtualVerletList");
  
    if (!connectionResort.connected()) {
      connectionResort.disconnect();
    }
  }
  
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VirtualVerletList::registerPython() {
    using namespace espresso::python;

    bool (VirtualVerletList::*pyExclude)(longint pid1, longint pid2)
          = &VirtualVerletList::exclude;


    class_<VirtualVerletList, shared_ptr<VirtualVerletList> >
      ("VirtualVerletList", init< shared_ptr<System>, real, shared_ptr<FixedTupleList>, bool >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VirtualVerletList::getBuilds, &VirtualVerletList::setBuilds)
      .def("totalSize", &VirtualVerletList::totalSize)
      .def("localSize", &VirtualVerletList::localSize)
      .def("getPair", &VirtualVerletList::getPair)
      .def("exclude", pyExclude)
      .def("rebuild", &VirtualVerletList::rebuild)
      .def("connect", &VirtualVerletList::connect)
      .def("disconnect", &VirtualVerletList::disconnect)
    
      .def("getVerletCutoff", &VirtualVerletList::getVerletCutoff)
      .def("setCellList", &VirtualVerletList::setCellList)
      .def("addTypeToVLMap", &VirtualVerletList::addTypeToVLMap)
      ;
  }

}
