/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"
#include "Collectives.hpp"

using namespace boost::mpi;
using namespace espressopp::python;
using namespace espressopp::esutil::Collectives;


DuplicateError::DuplicateError():
      std::runtime_error("item was found on more than one node")
{}

/**
   Reduction operator that evaluates to
   - NotHere iff all input values are NotHere
   - Duplicate iff more than one input value is not NotHere
   - the one value that is not NotHere otherwise
 */
struct UniqueReduce: public std::binary_function< int, int, int > {
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

 
int espressopp::esutil::Collectives::locateItem(bool here, int controller, communicator world) {
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
  return locateItem(here, controller, *mpiWorld);
}

namespace boost { namespace mpi {
  template<>
  struct is_commutative< UniqueReduce, int >: mpl::true_ { };
} } 

void espressopp::esutil::Collectives::registerPython() {
  def("esutil_Collectives_locateItem", pyLocateItem);
  scope().attr("esutil_Collectives_ResultNone") = None;
}
