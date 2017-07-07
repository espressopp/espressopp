/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
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
#ifndef _FIXEDTUPLELISTADRESS_HPP
#define _FIXEDTUPLELISTADRESS_HPP

#include "log4espp.hpp"
//#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espressopp {
    class FixedTupleListAdress: public TupleList  {
        protected:
            boost::signals2::connection sigOnTupleChanged, sigAfterRecv, sigBeforeSend;
            shared_ptr<storage::Storage> storage;
            typedef std::vector<longint> tuple;
            typedef boost::unordered_map<longint, tuple> GlobalTuples;
            GlobalTuples globalTuples;
            using TupleList::add;

        public:
            FixedTupleListAdress(shared_ptr<storage::Storage> _storage);
            ~FixedTupleListAdress();

            void add(longint pid) { tmppids.push_back(pid); } // add particle id (called from python)
            void addTs() { addT(tmppids); tmppids.clear(); } // add tuple (called from python)
            void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
            void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
            void onParticlesChanged();
            
            //int getNumPart(longint pid); // get number of particles in globalmap for given pid

            // This signals the AT particles to rebuild AT fixed pair, triple, quadruple bonds
            // Fixed{Pair|Triple}ListAdress connects to it
            boost::signals2::signal<void (std::vector<longint>&, class OutBuffer&)>
                beforeSendATParticles;
            boost::signals2::signal<void (ParticleList&, class InBuffer&)>
                afterRecvATParticles;

            static void registerPython();

        private:
            tuple tmppids;
            bool addT(tuple pids); // add tuple
            static LOG4ESPP_DECL_LOGGER(theLogger);
    };
}

#endif

