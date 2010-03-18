#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"
#include "boost/enable_shared_from_this.hpp"
#include "interaction/Interaction.hpp"

namespace espresso {

  namespace esutil { class RNG; }

  class System : public enable_shared_from_this<System> {

  public:
    System() : name("DEFAULT") {};
    System(std::string &_name) : name(_name) {};

    std::string name;
    shared_ptr< storage::Storage > storage;
    shared_ptr< bc::BC > bc;
    shared_ptr< esutil::RNG > rng;

    interaction::InteractionList shortRangeInteractions;

    real skin;  //<! skin used for VerletList

    shared_ptr<System> getShared() { 
      return shared_from_this();
    }

    static void registerPython();

    //!< keeps a pointer back to the Python object

  };
}
#endif
