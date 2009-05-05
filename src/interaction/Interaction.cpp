
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include "interaction/Interaction.hpp"

using namespace espresso::interaction;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Interaction::registerPython() {

  using namespace boost::python;

  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_<Interaction, boost::shared_ptr<Interaction>, boost::noncopyable >("interaction_Interaction", no_init)
  ;
}

