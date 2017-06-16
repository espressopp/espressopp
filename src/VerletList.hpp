/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)s
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
#ifndef _VERLETLIST_HPP
#define _VERLETLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"
#include "FixedPairList.hpp"
#include "FixedTripleList.hpp"
#include "FixedQuadrupleList.hpp"
#include "esutil/Timer.hpp"

namespace espressopp {
typedef boost::unordered_set<std::pair<longint, longint> > ExcludeList;

class DynamicExcludeList {
 public:
  DynamicExcludeList(shared_ptr<integrator::MDIntegrator> integrator);
  ~DynamicExcludeList();
  void exclude(longint pid1, longint pid2);
  void unexclude(longint pid1, longint pid2);
  void connect();
  void disconnect();
  shared_ptr<ExcludeList> getExList() { return exList; };
  python::list getList();
  int getSize() const { return exList->size(); }

  void observe_tuple(shared_ptr<FixedPairList> fpl);
  void observe_triple(shared_ptr<FixedTripleList> ftl);
  void observe_quadruple(shared_ptr<FixedQuadrupleList> fql);

  boost::signals2::signal<void ()> onListUpdated;
  boost::signals2::signal<void (longint, longint)> onPairExclude;
  boost::signals2::signal<void (longint, longint)> onPairUnexclude;
  
  static void registerPython();

 private:
  shared_ptr<integrator::MDIntegrator> integrator_;
  shared_ptr<System> system_;
  shared_ptr<ExcludeList> exList;
  // Helper lists.
  std::vector<longint> exList_add;
  std::vector<longint> exList_remove;

  bool is_dirty_;

  /**
   * Synchronize exclude list among all CPUs.
   */
  void updateList();

  boost::signals2::connection befIntP, runInit;
  static LOG4ESPP_DECL_LOGGER(theLogger);

};

/** Class that builds and stores verlet lists.

    ToDo: register at system for rebuild

*/

  class VerletList : public SystemAccess {

  public:

    /** Build a verlet list of all particle pairs in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list is built
	\param cut is the cutoff value for the 

    */

    VerletList(shared_ptr< System >, real cut, bool rebuildVL);
    VerletList(shared_ptr< System >, real cut, shared_ptr<DynamicExcludeList> dynamicExList, bool rebuildVL);

    ~VerletList();

    PairList& getPairs() { return vlPairs; }

    python::tuple getPair(int i);
    
    real getVerletCutoff(); // returns cutoff + skin

    void setVerletCutoff(real _cut); // set cutoff

    void connect();

    void disconnect();

    void rebuild();

    /** Get the total number of pairs for the Verlet list */
    int totalSize() const;

    //** Get the number of pairs for the local Verlet list */
    int localSize() const;

    /** Add pairs to exclusion list */
    bool exclude(longint pid1, longint pid2);

    /** Remove pairs from exclusion list. */
    bool unexclude(longint pid1, longint pid2);

    longint excludeListSize() const;

    /** Get the number of times the Verlet list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

    boost::signals2::signal<void (longint, longint)> onPairExclude;
    boost::signals2::signal<void (longint, longint)> onPairUnexclude;

  protected:

    void checkPair(Particle &pt1, Particle &pt2);
    PairList vlPairs;
    shared_ptr<ExcludeList> exList; // exclusion list
    shared_ptr<DynamicExcludeList> dynamicExcludeList;
    bool isDynamicExList;
    
    real cutsq;
    real cut;
    real cutVerlet;
    
    int builds;
    boost::signals2::connection connectionResort;

    /** timers */
    esutil::WallTimer wallTimer;
    real timeRebuild_;
    python::list getTimers() {
      python::list ret;
      ret.append(python::make_tuple("timeRebuild", timeRebuild_));
      return ret;
    }

    void resetTimers() {
      timeRebuild_ = 0.0;
      wallTimer.reset();
    }


    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
