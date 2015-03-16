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
#ifndef _VERLETLISTTRIPLE_HPP
#define _VERLETLISTTRIPLE_HPP

#include "python.hpp"
#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

namespace espressopp {

/** Class that builds and stores verlet lists for 3-body interactions.
    ToDo: register at system for rebuild
*/

  class VerletListTriple : public SystemAccess {

  public:

    /** Build a verlet list of all particle triples in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list triples are built
	\param cut is the cutoff value for the 

    */

    VerletListTriple(shared_ptr< System >, real cut, bool rebuildVL);

    ~VerletListTriple();

    TripleList& getTriples() { return vlTriples; }

    python::tuple getTriple(int i);

    real getVerletCutoff(); // returns cutoff + skin

    void connect();

    void disconnect();

    void rebuild();

    /** Get the total number of triples for the Verlet Triple list */
    int totalSize() const;

    //** Get the number of triples for the local Verlet list */
    int localSize() const;
    
    /** Add particle to exclusion list */
    bool exclude(longint pid);

    /** Get the number of times the Verlet Triple list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet Triple  list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

  protected:

    void checkTriple(Particle &pt1, Particle &pt2, Particle &pt3);
    TripleList vlTriples;
    
    boost::unordered_set< longint> exList; // exclusion list
    
    real cutsq;
    real cut;
    real cutVerlet;
    
    int builds;
    boost::signals2::connection connectionResort;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
