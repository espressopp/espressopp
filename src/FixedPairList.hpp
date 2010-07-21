// ESPP_CLASS
#ifndef _FIXEDPAIRLIST_HPP
#define _FIXEDPAIRLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espresso {
  class FixedPairList : public PairList {
  protected:
    boost::signals2::connection con1, con2, con3;
    shared_ptr< storage::Storage > storage;
    typedef boost::unordered_multimap< longint, longint > GlobalPairs;
    GlobalPairs globalPairs;

  public:
    using PairList::add;

    FixedPairList(shared_ptr< storage::Storage > _storage);
    ~FixedPairList();

    bool add(longint pid1, longint pid2);

    void beforeSendParticles(ParticleList& pl, 
			     mpi::packed_oarchive& ar);
    void afterRecvParticles(ParticleList& pl, 
			    mpi::packed_iarchive& ar);
    void onResortParticles();
    static void registerPython();
  private:
    static LOG4ESPP_DECL_LOGGER(theLogger);
  };
}

#endif

