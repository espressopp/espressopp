#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"
#include "Interaction.hpp"

namespace espresso {

  namespace esutil { class RNG; }

  class System {

  public:
    System() : name("DEFAULT") {};
    System(std::string &_name) : name(_name) {};

    std::string name;
    shared_ptr< Storage > storage;
    shared_ptr< BC > bc;
    shared_ptr< esutil::RNG > rng;

    interaction::InteractionList shortRangeInteractions;

    double skin;  //<! skin used for VerletList

    static void registerPython();

  };
}
#endif
