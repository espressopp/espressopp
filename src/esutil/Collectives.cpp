#include "Collectives.hpp"
#include <boost/python.hpp>

using namespace boost::mpi;
using namespace boost::python;
using namespace espresso::esutil::Collectives;

/**
   Reduction operator that evaluates to
   - NotHere iff all input values are NotHere
   - Duplicate iff more than one input value is not NotHere
   - the one value that is not NotHere otherwise
 */
struct UniqueReduce: public std::binary_function<int, int, int> {
  static const int NotHere   = -1;
  static const int Duplicate = -2;

  int operator() (int x, int y) {
    if (x == NotHere) {
      return y;
    }
    else if (y == NotHere) {
      return x;
    }
    else {
      return Duplicate;
    }
  }
};

DuplicateError::DuplicateError():
      std::runtime_error("item was found on more than one node")
{}
 
namespace boost { namespace mpi {
  template<>
  struct is_commutative<UniqueReduce, int>: mpl::true_ { };
} } 

int espresso::esutil::Collectives::locateItem(bool here, int controller, communicator world) {
  int node = here ? world.rank() : UniqueReduce::NotHere;

  if (world.rank() != controller) {
    reduce(world, node, UniqueReduce(), controller);
    return None;
  }
  else {
    int owner;
    reduce(world, node, owner, UniqueReduce(), controller);
    if (owner == UniqueReduce::Duplicate) {
      throw DuplicateError();
    }
    return owner;
  }
}

static int pyLocateItem(bool here, int controller) {
  return locateItem(here, controller, communicator());
}

void espresso::esutil::Collectives::registerPython() {
  
  def("esutil_Collectives_locateItem", pyLocateItem);
  scope().attr("esutil_Collectives_ResultNone") = None;

}
