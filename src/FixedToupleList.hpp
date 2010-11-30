// ESPP_CLASS
#ifndef __ESPRESSOPP_FIXEDPAIRLIST_HPP
#define __ESPRESSOPP_FIXEDPAIRLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espresso {
  //template<typename ToupleType = std::vector<longint> > // we need to add creator policy

  /**
   * Start implementing virtual sites here, will be generalized later on
   */class FixedToupleList  {
  protected:
    boost::signals2::connection con1, con2, con3;
    shared_ptr< storage::Storage > storage;
    //typedef boost::unordered_multimap< longint, longint > GlobalPairs;
    //GlobalPairs globalPairs;

    //using PairList::add;

  public:
    FixedToupleList(shared_ptr< storage::Storage > _storage);
    ~FixedToupleList();

    //bool add(longint pid1, longint pid2);

    void beforeSendParticles(ParticleList& pl, 
			     class OutBuffer& buf);
    void afterRecvParticles(ParticleList& pl, 
			    class InBuffer& buf);
    void onParticlesChanged();

    static void registerPython();
  private:
    static LOG4ESPP_DECL_LOGGER(theLogger);
  };
}

#endif

