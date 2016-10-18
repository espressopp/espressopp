/*
  Copyright (C) 2012,2013,2014,2015,2016
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
#ifndef _ANALYSIS_RDFATOMISTIC_HPP
#define _ANALYSIS_RDFATOMISTIC_HPP

#include "types.hpp"
#include "Observable.hpp"

#include "python.hpp"

namespace espressopp {
  namespace analysis {
    /** Class to compute the radial distribution function of the system. */
    class RDFatomistic : public Observable {
    public:
      RDFatomistic(shared_ptr< System > system, int type1, int type2, bool _spanbased = true, real _span = 1.0) : Observable(system), target1(type1), target2(type2), spanbased(_spanbased), span(_span) {}
      ~RDFatomistic() {}
      virtual real compute() const;
      virtual python::list computeArray(int) const;

      static void registerPython();

      class data {
        public:
        Real3D pos;
        int type;
        int molecule;
        real resolution;
        friend class boost::serialization::access;
        template<class Archive> void serialize(Archive & ar, const unsigned int version) {
          ar & pos;
          ar & type;
          ar & molecule;
          ar & resolution;
        }
      };

      int target1;
      int target2;
      real span;
      bool spanbased;

    };
  }
}

#endif
