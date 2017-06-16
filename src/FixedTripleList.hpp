/*
  Copyright (C) 2012,2013,2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
  
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
#ifndef _FIXEDTRIPLELIST_HPP
#define _FIXEDTRIPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espressopp {
  class FixedTripleList : public TripleList {
      protected:
		boost::signals2::connection sigAfterRecv, sigOnParticleChanged, sigBeforeSend;
		shared_ptr<storage::Storage> storage;
		typedef boost::unordered_multimap <longint,std::pair <longint, longint> > GlobalTriples;
		GlobalTriples globalTriples;
		using TripleList::add;

      //FixedListComm<FixedTripleList, 3> _comm;

	  public:
	    FixedTripleList() { }
		FixedTripleList(shared_ptr<storage::Storage> _storage);
		virtual ~FixedTripleList();
		//bool add(pvec pids) { _comm.add(pids); }
		/** Add the given particle triple to the list on this processor if the
		particle with the lower id belongs to this processor.  Note that
		this routine does not check whether the triple is inserted on
		another processor as well.
		\return whether the triple was inserted on this processor.
		*/
		virtual bool add(longint pid1, longint pid2, longint pid3);
	    // Non blocking version of 'add' method.
	    virtual bool iadd(longint pid1, longint pid2, longint pid3);

        /**
         * Removes a triplet from the list.
         * @param pid1 particle id
         * @param pid2 particle id
         * @param pid3 particle id
         * @param no_signal if true, onTupleRemoved signal will not be thrown
         * @return true if the triplet is removed
         */
		virtual bool remove(longint pid1, longint pid2, longint pid3, bool no_signal);
		virtual void clearAndRemove();
		virtual bool removeByBond(longint pid1, longint pid2);

		virtual void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
		virtual void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
		virtual void onParticlesChanged();
		virtual void updateParticlesStorage();

		virtual std::vector<longint> getTripleList();
		virtual python::list getTriples();
		virtual python::list getAllTriples();
	    /** Get the number of triples in the GlobalTriples list */
	    virtual int size() { return globalTriples.size(); }
	    virtual int totalSize();

	    boost::signals2::signal<void (longint, longint, longint)> onTupleAdded;
	    boost::signals2::signal<void (longint, longint, longint)> onTupleRemoved;

		static void registerPython();

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);
  };
}

#endif
