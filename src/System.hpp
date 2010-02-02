#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  namespace interaction { class Interaction; }

  typedef std::vector< shared_ptr<interaction::Interaction> > InteractionList;

  class System {

  public:
    shared_ptr< class Storage > storage;
    shared_ptr< class BC > bc;

    InteractionList shortRangeInteractions;

    esutil::RNG rng; //<! common random number generator

    double skin;  //<! skin used for VerletList

  };
}
#endif
