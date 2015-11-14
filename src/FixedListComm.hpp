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
#ifndef _FIXEDLISTCOMM_HPP
#define _FIXEDLISTCOMM_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
//#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espressopp {
    class FixedListComm: public PairList, public TripleList, public QuadrupleList {
        protected:
        boost::signals2::connection con1, con2, con3;
        shared_ptr<storage::Storage> storage;
        typedef std::vector<longint> pvec;
        typedef boost::unordered_multimap<longint, pvec> GlobalList;
        GlobalList globalLists;
        //using PairList::add;

        // if one wants to get rid of the virtual,
        // move this to a template
        //virtual void add(std::vector<Particle*> tmp)=0;

        public:
        FixedListComm(shared_ptr<storage::Storage> _storage);
        ~FixedListComm();
        bool add(pvec pids);
        void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
        void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
        void onParticlesChanged();

        private:
        static LOG4ESPP_DECL_LOGGER(theLogger);
    }; // class
} //ns espressopp

#endif
