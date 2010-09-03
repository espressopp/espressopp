// ESPP_CLASS
#ifndef _FIXEDTRIPLELIST_HPP
#define _FIXEDTRIPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espresso {
  class FixedTripleList : public TripleList {
  protected:
    boost::signals2::connection con1, con2, con3;
    shared_ptr< storage::Storage > storage;
    typedef boost::unordered_multimap< longint, std::pair < longint, longint > > GlobalTriples;
    GlobalTriples globalTriples;

    using TripleList::add;

  public:
    FixedTripleList(shared_ptr< storage::Storage > _storage);
    ~FixedTripleList();

    /** Add the given particle triple to the list on this processor if the
	particle with the lower id belongs to this processor.  Note that
	this routine does not check whether the triple is inserted on
	another processor as well.  
	
	\return whether the triple was inserted on this processor.
    */
    bool add(longint pid1, longint pid2, longint pid3);

    void beforeSendParticles(ParticleList& pl, 
			     mpi::packed_oarchive& ar);
    void afterRecvParticles(ParticleList& pl, 
			    mpi::packed_iarchive& ar);
    void onParticlesChanged();
    static void registerPython();
  private:
    static LOG4ESPP_DECL_LOGGER(theLogger);
  };
}

#endif
