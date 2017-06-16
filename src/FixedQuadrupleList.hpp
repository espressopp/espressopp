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
#ifndef _FIXEDQUADRUPLELIST_HPP
#define _FIXEDQUADRUPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Triple.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espressopp {
  class FixedQuadrupleList : public QuadrupleList {
  protected:
    boost::signals2::connection sigBeforeSend, sigAfterRecv, sigOnParticlesChanged;
    shared_ptr< storage::Storage > storage;
    typedef boost::unordered_multimap< longint,
            Triple < longint, longint, longint > > GlobalQuadruples;
    GlobalQuadruples globalQuadruples;
    using QuadrupleList::add;

  public:
    FixedQuadrupleList() { }
    FixedQuadrupleList(shared_ptr< storage::Storage > _storage);
    ~FixedQuadrupleList();

    /** Add the given particle quadruple to the list on this processor if the
	particle with the lower id belongs to this processor.  Note that
	this routine does not check whether the quadruple is inserted on
	another processor as well.  
	
	\return whether the quadruple was inserted on this processor.
    */
    virtual bool add(longint pid1, longint pid2, longint pid3, longint pid4);
    /// Non-blocking add method.
    virtual bool iadd(longint pid1, longint pid2, longint pid3, longint pid4);

    virtual bool remove(longint pid1, longint pid2, longint pid3, longint pid4, bool no_signal);
    virtual void clearAndRemove() {
      this->clear();
      globalQuadruples.clear();
      sigBeforeSend.disconnect();
      sigAfterRecv.disconnect();
    }
    /**
     * Remove quadruplets that are build around bond pid1-pid2
     * @param pid1 particle id
     * @param pid2 particle id
     * @return true if the quadruplets were removed
     */
    virtual bool removeByBond(longint pid1, longint pid2);

    virtual void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
    virtual void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
    virtual void onParticlesChanged();
    virtual void updateParticlesStorage();

    virtual std::vector<longint> getQuadrupleList();
    virtual python::list getQuadruples();

    /** Get the number of quadruples in the GlobalQuadruples list */
    virtual int size() {
    	return globalQuadruples.size();
    }

    virtual int totalSize();

    boost::signals2::signal<void (longint, longint, longint, longint)> onTupleAdded;
    boost::signals2::signal<void (longint, longint, longint, longint)> onTupleRemoved;

    static void registerPython();

  private:
    static LOG4ESPP_DECL_LOGGER(theLogger);
    python::list getAllQuadruples();
  };
}

#endif
