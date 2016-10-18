/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
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
#ifndef _ANALYSIS_ConfigurationsExt_HPP
#define _ANALYSIS_ConfigurationsExt_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "ConfigurationExt.hpp"

namespace espressopp {
  namespace analysis {

    /** Class that stores particle positions for later analysis.

        Important: this class can also be used if the number of
        particles changes between different snapshots.
    */

    typedef std::vector<ConfigurationExtPtr> ConfigurationExtList;

    class ConfigurationsExt : public SystemAccess {

    public:

      /** Constructor, allow for unlimited snapshots. */

      ConfigurationsExt(shared_ptr<System> system) : SystemAccess (system)
      { maxConfigs = 0; }

      /** set number of maximal snapshots. */

      void setCapacity(int max);

      /** get number of maximal snapshots. */

      int getCapacity();

      /** get number of available snapshots. */

      int getSize();

      ~ConfigurationsExt() {}

      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}

      /** Take a snapshot of all current particle positions. */

      void gather();

      ConfigurationExtPtr get(int stackpos);

      ConfigurationExtPtr back();

      ConfigurationExtList all();

      void clear() { configurationsExt.clear(); }

      static void registerPython();

    protected:

      static LOG4ESPP_DECL_LOGGER(logger);

    private:

      void pushConfig(ConfigurationExtPtr config);

      ConfigurationExtList configurationsExt;

      int maxConfigs;

      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is unfolded
    };
  }
}

#endif
