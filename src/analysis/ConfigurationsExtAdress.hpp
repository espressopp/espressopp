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
#ifndef _ANALYSIS_ConfigurationsExtAdress_HPP
#define _ANALYSIS_ConfigurationsExtAdress_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "ConfigurationExt.hpp"
#include "FixedTupleListAdress.hpp"

namespace espressopp {
  namespace analysis {

    /** Class that stores atomistic particle positions for later analysis.

        Important: this class can also be used if the number of
        particles changes between different snapshots.
    */

    typedef std::vector<ConfigurationExtPtr> ConfigurationExtList;

    class ConfigurationsExtAdress : public SystemAccess {

    public:

      /** Constructor, allow for unlimited snapshots. */

      ConfigurationsExtAdress(shared_ptr<System> system,shared_ptr<FixedTupleListAdress> _fixedTupleList) : SystemAccess (system),fixedTupleList(_fixedTupleList)
      { maxConfigs = 0; }

      /** set number of maximal snapshots. */

      void setCapacity(int max);

      /** get number of maximal snapshots. */

      int getCapacity();

      /** get number of available snapshots. */

      int getSize();

      ~ConfigurationsExtAdress() {}
      
      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}

      /** Gake a snapshot of all current particle positions. */

      void gather();

      ConfigurationExtPtr get(int stackpos);

      ConfigurationExtPtr back();

      ConfigurationExtList all();

      void clear() { configurationsExtAdress.clear(); }

      void
      setFixedTupleList(shared_ptr<FixedTupleListAdress> _fixedtupleList) {
          fixedTupleList = _fixedtupleList;
      }

      static void registerPython();
    
    protected:

      static LOG4ESPP_DECL_LOGGER(logger);
      shared_ptr<FixedTupleListAdress> fixedTupleList;

    private:

      void pushConfig(ConfigurationExtPtr config);
 
      ConfigurationExtList configurationsExtAdress;

      int maxConfigs;
      
      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is unfolded
    };
  }
}

#endif
