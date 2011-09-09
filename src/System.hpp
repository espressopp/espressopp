// ESPP_CLASS
#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"
#include "mpi.hpp"
#include "python.hpp"
#include "boost/enable_shared_from_this.hpp"
#include "interaction/Interaction.hpp"

namespace espresso {

  namespace esutil { class RNG; }

  class System : public enable_shared_from_this< System > {

  public:

    System();
    System(python::object _pyobj);

    shared_ptr< mpi::communicator > comm;

    shared_ptr< storage::Storage > storage;
    shared_ptr< bc::BC > bc;
    shared_ptr< esutil::RNG > rng;

    interaction::InteractionList shortRangeInteractions;

    real skin;  //<! skin used for VerletList

    bool CommunicatorIsInitialized;

    shared_ptr< System > getShared() { 
      return shared_from_this();
    }

    void addInteraction(shared_ptr< interaction::Interaction > ia);
    shared_ptr< interaction::Interaction > getInteraction(int i);
    int getNumberOfInteractions();
    static void registerPython();

  };
}
#endif
