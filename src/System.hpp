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

// ESPP_CLASS
#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "python.hpp"
#include "mpi.hpp"
#include "boost/enable_shared_from_this.hpp"
#include "interaction/Interaction.hpp"
#include "types.hpp"

namespace espressopp {

  namespace esutil { 
    class RNG;
  }

  class System : public enable_shared_from_this< System > {
  
  private:
    real skin;  //<! skin used for VerletList
    
  public:

    System();
    System(python::object _pyobj);

    shared_ptr< mpi::communicator > comm;

    shared_ptr< storage::Storage > storage;
    shared_ptr< bc::BC > bc;
    shared_ptr< esutil::RNG > rng;

    interaction::InteractionList shortRangeInteractions;

    real maxCutoff;     // maximal cutoff over all of the interactions

    bool CommunicatorIsInitialized;

    shared_ptr< System > getShared() { 
      return shared_from_this();
    }
    
    void setSkin(real);
    real getSkin();
    
    void scaleVolume(real s, bool particleCoordinates);
    void scaleVolume(Real3D s, bool particleCoordinates);
    void scaleVolume3D(Real3D s);
    void setTrace(bool flag);
    void addInteraction(shared_ptr< interaction::Interaction > ia);
    void removeInteraction(int i);
    shared_ptr< interaction::Interaction > getInteraction(int i);
    int getNumberOfInteractions();
    static void registerPython();

  };
}
#endif
