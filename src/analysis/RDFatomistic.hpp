/*
  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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
        RDFatomistic(shared_ptr< System > system, int type1, int type2, real _span = 1.0, bool _spanbased = true) : Observable(system), target1(type1), target2(type2), span(_span), spanbased(_spanbased) {}
        ~RDFatomistic() {}
        virtual real compute() const;
        virtual python::list computeArray(int) const;
        virtual python::list computeArrayPathIntegral(int) const;

        static void registerPython();

        class data {
          public:
            Real3D pos;
            size_t type;
            size_t molecule;
            real resolution;
            friend class boost::serialization::access;
            template<class Archive> void serialize(Archive & ar, const unsigned int version) {
              ar & pos;
              ar & type;
              ar & molecule;
              ar & resolution;
            }
        };

        class dataPathIntegral {
          public:
            Real3D pos;
            size_t type;
            size_t molecule;
            size_t pib;
            friend class boost::serialization::access;
            template<class Archive> void serialize(Archive & ar, const unsigned int version) {
              ar & pos;
              ar & type;
              ar & molecule;
              ar & pib;
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
