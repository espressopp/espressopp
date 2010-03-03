#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"

namespace espresso {

  namespace interaction { class Interaction; }
  namespace esutil { class RNG; }

  typedef std::vector< shared_ptr<interaction::Interaction> > InteractionList;

  class System {

  public:
    System() : name("DEFAULT") {};
    System(std::string _name) : name(_name) {};

    std::string name;
    shared_ptr< class Storage > storage;
    shared_ptr< class BC > bc;
    shared_ptr< esutil::RNG > rng;

    InteractionList shortRangeInteractions;

    double skin;  //<! skin used for VerletList

    static void registerPython();

  };
}
#endif
