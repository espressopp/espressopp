// ESPP_CLASS
#ifndef _FIXEDQUADRUPLELIST_HPP
#define _FIXEDQUADRUPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Triple.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espresso {
  class FixedQuadrupleList : public QuadrupleList {
  protected:
    boost::signals2::connection con1, con2, con3;
    shared_ptr< storage::Storage > storage;
    typedef boost::unordered_multimap< longint,
            Triple < longint, longint, longint > > GlobalQuadruples;
    GlobalQuadruples globalQuadruples;

    using QuadrupleList::add;

  public:
    FixedQuadrupleList(shared_ptr< storage::Storage > _storage);
    ~FixedQuadrupleList();

    /** Add the given particle quadruple to the list on this processor if the
	particle with the lower id belongs to this processor.  Note that
	this routine does not check whether the quadruple is inserted on
	another processor as well.  
	
	\return whether the quadruple was inserted on this processor.
    */
    bool add(longint pid1, longint pid2, longint pid3, longint pid4);

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
