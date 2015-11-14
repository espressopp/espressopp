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
#ifndef _FIXEDQUADRUPLEANGLELIST_HPP
#define _FIXEDQUADRUPLEANGLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Triple.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
//#include <boost/unordered_map.hpp>
#include <map>
#include <boost/signals2.hpp>

namespace espressopp {
  class FixedQuadrupleAngleList : public QuadrupleList{
  protected:
    boost::signals2::connection con1, con2, con3;
    typedef std::multimap< longint,
            std::pair<Triple < longint, longint, longint >, real> > QuadruplesAngles;
    shared_ptr <storage::Storage> storage;
    QuadruplesAngles quadruplesAngles;
    using QuadrupleList::add;

  public:
    FixedQuadrupleAngleList(shared_ptr <storage::Storage> _storage);
    ~FixedQuadrupleAngleList();

    /** Add the given particle quadruple to the list on this processor if the
	particle with the lower id belongs to this processor.  Note that
	this routine does not check whether the quadruple is inserted on
	another processor as well.  
     * 
     * it contains the information about the initial angle for each quadruple.
     * 
     * This is probably temporary fast solution for the dihedral potential where 
     * all the quadruples have different phi0 
     * 
	
	\return whether the quadruple was inserted on this processor.
    */
    bool add(longint pid1, longint pid2, longint pid3, longint pid4);
    void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
    void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
    void onParticlesChanged();

    python::list getQuadruples();
	// get triples and corresponding angles
	python::list getQuadruplesAngles();

    // get angle value for current triple
    real getAngle(int, int, int, int);
    /** Get the number of quadruples in the GlobalQuadruples list */
    int size() {return quadruplesAngles.size();}

    static void registerPython();

  private:
    static LOG4ESPP_DECL_LOGGER(theLogger);
  };
}

#endif
