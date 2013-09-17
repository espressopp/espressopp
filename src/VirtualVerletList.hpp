// ESPP_CLASS
#ifndef _VirtualVerletList_HPP
#define _VirtualVerletList_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"
#include "Cell.hpp"
#include "VerletList.hpp"
#include "esutil/Array2D.hpp"
#include "FixedTupleList.hpp"

namespace espresso {

/** Class that builds and stores verlet lists.

    ToDo: register at system for rebuild

*/

  class VirtualVerletList : public SystemAccess {

  public:

    /** Build a verlet list of all particle pairs in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list is built
	\param cut is the cutoff value for the 

    */

    VirtualVerletList(shared_ptr< System >, real cut, shared_ptr<FixedTupleList> ftpl, bool rebuildVL);

    ~VirtualVerletList();

    PairList& getPairs() { return vlPairs; }

    python::tuple getPair(int i);
    
    real getVerletCutoff(); // returns cutoff + skin

    void connect();

    void disconnect();

    void rebuild();

    /** Get the total number of pairs for the Verlet list */
    int totalSize() const;

    //** Get the number of pairs for the local Verlet list */
    int localSize() const;

    /** Add pairs to exclusion list */
    bool exclude(longint pid1, longint pid2);

    /** Get the number of times the Verlet list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

    void addTypeToVLMap(int typeI, int typeJ, shared_ptr<VerletList> vl){
    	vlArray.at(typeI, typeJ)= vl;
    }

    void setCellList(shared_ptr<CellList> _cellList){
    	std::cout<<"cell list set " << _cellList->size() << std::endl;
    	cellList=_cellList;
    };


  protected:

    void checkPair(Particle &pt1, Particle &pt2);
    PairList vlPairs;
    boost::unordered_set<std::pair<longint, longint> > exList; // exclusion list
    
    real cutsq;
    real cut;
    real cutVerlet;
    
    int builds;
    boost::signals2::connection connectionResort;

    shared_ptr<CellList> cellList;

    typedef std::map<int, shared_ptr<VerletList> > TypeToVLMap;
    TypeToVLMap typemap;

    esutil::Array2D<shared_ptr<VerletList>, esutil::enlarge> vlArray;

    shared_ptr <FixedTupleList> ftpl;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
