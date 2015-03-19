/*
  Copyright (C) 2014
      Pierre de Buyl
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

#include "python.hpp"
#include "AssociationReaction.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"

namespace espressopp {

  namespace integrator {

    //using namespace espressopp::iterator;
    using namespace storage;

    LOG4ESPP_LOGGER(AssociationReaction::theLogger, "AssociationReaction");

    AssociationReaction::AssociationReaction(
					     shared_ptr<System> system,
					     shared_ptr<VerletList> _verletList,
					     shared_ptr<FixedPairList> _fpl,
					     shared_ptr<DomainDecomposition> _domdec)
      :Extension(system), verletList(_verletList), fpl(_fpl), domdec(_domdec) {

      type = Extension::Reaction;

      rate  = 0.0;
      cutoff = 0.0;
      typeA = 0;
      typeB = 0;
      deltaA = 0;
      deltaB = 0;
      stateAMin = 0;

      current_cutoff = verletList->getVerletCutoff() - system->getSkin();
      current_cutoff_sqr = current_cutoff*current_cutoff;
      
      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "AssociationReaction constructed");
    }

    void AssociationReaction::setRate(real _rate)
    {
      rate = _rate;
    }

    real AssociationReaction::getRate()
    {
      return rate;
    }

    void AssociationReaction::setCutoff(real _cutoff)
    {
      cutoff = _cutoff;
      cutoff_sqr = _cutoff*_cutoff;
    }

    real AssociationReaction::getCutoff()
    {
      return cutoff;
    }

    void AssociationReaction::setTypeA(size_t _typeA)
    {
      typeA = _typeA;
    }

    size_t AssociationReaction::getTypeA()
    {
      return typeA;
    }

    void AssociationReaction::setTypeB(size_t _typeB)
    {
      typeB = _typeB;
    }

    size_t AssociationReaction::getTypeB()
    {
      return typeB;
    }

    void AssociationReaction::setDeltaA(int _deltaA)
    {
      deltaA = _deltaA;
    }

    int AssociationReaction::getDeltaA()
    {
      return deltaA;
    }

    void AssociationReaction::setDeltaB(int _deltaB)
    {
      deltaB = _deltaB;
    }

    int AssociationReaction::getDeltaB()
    {
      return deltaB;
    }

    void AssociationReaction::setStateAMin(int _stateAMin)
    {
      stateAMin = _stateAMin;
    }

    int AssociationReaction::getStateAMin()
    {
      return stateAMin;
    }

    void AssociationReaction::setInterval(int _interval)
    {
      interval = _interval;
    }

    int AssociationReaction::getInterval()
    {
      return interval;
    }

    AssociationReaction::~AssociationReaction() {
      disconnect();
    }

    void AssociationReaction::disconnect() {

      _initialize.disconnect();
      _react.disconnect();
    }

    void AssociationReaction::connect() {

      // connect to initialization inside run()
      _initialize = integrator->runInit.connect(
						boost::bind(&AssociationReaction::initialize, this));

      _react = integrator->aftIntV.connect(
					   boost::bind(&AssociationReaction::react, this));
    }

    /** Performs all steps of the reactive scheme.
     */
    void AssociationReaction::react() {
      if (integrator->getStep() % interval != 0) return;

      System& system = getSystemRef();

      LOG4ESPP_INFO(theLogger, "Perform AssociationReaction");

      dt = integrator->getTimeStep();

      Alist.clear();
      // loop over VL pairs
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
	Particle &p1 = *it->first;
	Particle &p2 = *it->second;
	// If criteria for reaction match, add the indices to Alist
	reactPair(p1, p2);
      }
      sendMultiMap(Alist);
      // Here, reduce number of partners to each A to 1
      // Also, keep only non-ghost A
      uniqueA(Alist);
      sendMultiMap(Alist);
      // Here, reduce number of partners to each B to 1
      // Also, keep only non-ghost B
      uniqueB(Alist, Blist);
      sendMultiMap(Blist);
      // Use Blist to apply the reaction.
      applyAR();
    }

    /** For a given pair of particles, check if they meet the condition the
	reactive scheme. If it is so, the (A,B) pair is added to Alist.
    */
    void AssociationReaction::reactPair(Particle& p1, Particle& p2) {
      Real3D r = p1.position() - p2.position();
      real dist2 = r.sqr();
      bool found;
      found=false;
      if ((dist2 < cutoff_sqr) && ((*rng)() < rate*dt*interval)){
	if ((p1.type()==typeA) && (p2.type()==typeB) && (p1.state() >= stateAMin) && (p2.state()==0)) {
	  Alist.insert(std::make_pair(p1.id(), p2.id()));
	}
	else if ((p2.type()==typeA) && (p1.type()==typeB) && (p2.state() >= stateAMin) && (p1.state()==0)) {
	  Alist.insert(std::make_pair(p2.id(), p1.id()));
	}
      }
    }

    void AssociationReaction::initialize() {

      LOG4ESPP_INFO(theLogger, "init AssociationReaction");

    }

    /** Performs two-way parallel communication to consolidate mm between
	neighbours. The parallel scheme is taken from
	DomainDecomposition::doGhostCommunication
    */
    void AssociationReaction::sendMultiMap(boost::unordered_multimap<longint, longint> &mm) {

      LOG4ESPP_INFO(theLogger, "Entering sendMultiMap");

      InBuffer inBuffer0(*getSystem()->comm);
      InBuffer inBuffer1(*getSystem()->comm);
      OutBuffer outBuffer(*getSystem()->comm);
      System& system = getSystemRef();
      const NodeGrid& nodeGrid = domdec->getNodeGrid();

      /* direction loop: x, y, z.
	 Here we could in principle build in a one sided ghost
	 communication, simply by taking the lr loop only over one
	 value. */
      for (int coord = 0; coord < 3; ++coord) {
	/* inverted processing order for ghost force communication,
	   since the corner ghosts have to be collected via several
	   nodes. We now add back the corner ghost forces first again
	   to ghost forces, which only eventually go back to the real
	   particle.
	*/

	real curCoordBoxL = system.bc->getBoxL()[coord];

	outBuffer.reset();
	// fill outBuffer from mm
	int tmp = mm.size();
	outBuffer.write(tmp);
	for (boost::unordered_multimap<longint, longint>::iterator it=mm.begin(); it!=mm.end(); it++) {
	  tmp = it->first;
	  outBuffer.write(tmp);
	  tmp = it->second;
	  outBuffer.write(tmp);
	}
    
	// lr loop: left right
	for (int lr = 0; lr < 2; ++lr) {
	  int dir         = 2 * coord + lr;
	  int oppositeDir = 2 * coord + (1 - lr);
	  int dirSize = nodeGrid.getGridSize(coord);
	  // Avoids double communication for size 2 directions.
	  if ( (dirSize==2) && (lr==1) ) continue;

	  if (dirSize == 1) {
	    LOG4ESPP_DEBUG(theLogger, "no communication");
	  }
	  else {
	    // prepare send and receive buffers
	    longint receiver, sender;
	    receiver = nodeGrid.getNodeNeighborIndex(dir);
	    sender = nodeGrid.getNodeNeighborIndex(oppositeDir);

	    // exchange particles, odd-even rule
	    if (nodeGrid.getNodePosition(coord) % 2 == 0) {
	      outBuffer.send(receiver, AR_COMM_TAG);
	      if (lr==0) {
		inBuffer0.recv(sender, AR_COMM_TAG);
	      } else {
		inBuffer1.recv(sender, AR_COMM_TAG);
	      }
	    } else {
	      if (lr==0) {
		inBuffer0.recv(sender, AR_COMM_TAG);
	      } else {
		inBuffer1.recv(sender, AR_COMM_TAG);
	      }
	      outBuffer.send(receiver, AR_COMM_TAG);
	    }
	  }
	}
	LOG4ESPP_DEBUG(theLogger, "Entering unpack");
	// unpack received data
	// add content of inBuffer to mm
	int lengthA, Aidx, Bidx;
	for (int lr = 0; lr < 2; ++lr) {
	  int dir         = 2 * coord + lr;
	  int oppositeDir = 2 * coord + (1 - lr);
	  int dirSize = nodeGrid.getGridSize(coord);
	  if (dirSize == 1) {
	  } else {
	    // Avoids double communication for size 2 directions.
	    if ( (dirSize==2) && (lr==1) ) continue;
	    if (lr==0) {
	      inBuffer0.read(lengthA);
	    } else {
	      inBuffer1.read(lengthA);
	    }
	    for (longint i=0; i<lengthA; i++) {
	      if (lr==0) {
		inBuffer0.read(Aidx);
		inBuffer0.read(Bidx);
	      } else {
		inBuffer1.read(Aidx);
		inBuffer1.read(Bidx);
	      }
	      mm.insert(std::make_pair(Aidx, Bidx));
	    }
	  }
	}
	LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
      }

      LOG4ESPP_INFO(theLogger, "Leaving sendMultiMap");

    }

    /** Given a multimap mm with several pairs (id1,id2), keep only one pair for
	each id1 and return it in place. In addition, only pairs for which
	id1 is local are kept.
    */
    void AssociationReaction::uniqueA(boost::unordered_multimap<longint, longint> &mm) {

      // Collect indices
      boost::unordered_set<longint> idxSet;
      boost::unordered_multimap<longint, longint> idxList;
      longint idx1, idx2;
      System& system = getSystemRef();
      boost::unordered_multimap<longint, longint> uniqueList;

      idxSet.clear();
      idxList.clear();
      // Create idxList, containing (idx1, idx2) pairs
      // Create idxSet, containing each idx1 only once
      for (boost::unordered_multimap<longint, longint>::iterator it=mm.begin();
	   it!=mm.end(); it++) {
	idx1 = it->first;
	idx2 = it->second;
	idxList.insert(std::make_pair(idx1, idx2));
	idxSet.insert(idx1);
      }

      uniqueList.clear();
      // For each active idx1, pick a partner
      if (idxSet.size()>0) {
	for (boost::unordered_set<longint>::iterator it=idxSet.begin(); it!=idxSet.end(); it++) {
	  idx1 = *it;
	  Particle* p = system.storage->lookupLocalParticle(idx1);
	  if (p==NULL) continue;
	  if (p->ghost()) continue;
	  int size = idxList.count(idx1);
	  if (size>0) {
	    int pick = (*rng)(size);
	    std::pair<boost::unordered_multimap<longint, longint>::iterator,boost::unordered_multimap<longint, longint>::iterator> candidates;
	    candidates = idxList.equal_range(idx1);
	    int i=0;
	    for (boost::unordered_multimap<longint, longint>::iterator jt=candidates.first; jt!=candidates.second; jt++, i++) {
	      if (i==pick) {
		uniqueList.insert(std::make_pair(jt->first, jt->second));
		break;
	      }
	    }
	  }
	}
      }
      mm = uniqueList;
      uniqueList.clear();
    }

    /** Given a multimap mm with several pairs (id1,id2), keep only one pair for
	each id2 and return it in place. In addition, only pairs for which
	id2 is local are kept.
    */
    void AssociationReaction::uniqueB(boost::unordered_multimap<longint, longint> &mm, boost::unordered_multimap<longint, longint> &nn) {

      // Collect indices
      boost::unordered_set<longint> idxSet;
      boost::unordered_multimap<longint, longint> idxList;
      longint idx1, idx2;
      System& system = getSystemRef();
      boost::unordered_multimap<longint, longint> uniqueList;

      idxSet.clear();
      idxList.clear();
      // Create idxList, containing (idx2, idx1) pairs
      // Create idxSet, containing each idx2 only once
      for (boost::unordered_multimap<longint, longint>::iterator it=mm.begin();
	   it!=mm.end(); it++) {
	idx1 = it->first;
	idx2 = it->second;
	idxList.insert(std::make_pair(idx2, idx1));
	idxSet.insert(idx2);
      }

      uniqueList.clear();
      if (idxSet.size()>0) {
	// For each active idx1, pick a partner
	for (boost::unordered_set<longint>::iterator it=idxSet.begin(); it!=idxSet.end(); it++) {
	  idx2 = *it;
	  Particle* p = system.storage->lookupLocalParticle(idx2);
	  if (p==NULL) continue;
	  if (p->ghost()) continue;
	  int size = idxList.count(idx2);
	  if (size>0) {
	    int pick = (*rng)(size);
	    std::pair<boost::unordered_multimap<longint, longint>::iterator,boost::unordered_multimap<longint, longint>::iterator> candidates;
	    candidates = idxList.equal_range(idx2);
	    int i=0;
	    for (boost::unordered_multimap<longint, longint>::iterator jt=candidates.first; jt!=candidates.second; jt++, i++) {
	      if (i==pick) {
		uniqueList.insert(std::make_pair(jt->first, jt->second));
		break;
	      }
	    }
	  }
	}
      }
      nn = uniqueList;
      uniqueList.clear();
    }

    /** Use the (A,B) list "partners" to add bonds and change the state of the
	particles accordingly.
    */
    void AssociationReaction::applyAR() {
      longint A, B;
      System& system = getSystemRef();

      LOG4ESPP_INFO(theLogger, "Entering applyAR");

      for (boost::unordered_multimap<longint, longint>::iterator it=Blist.begin(); it!=Blist.end(); it++) {
	A = it->second;
	B = it->first;
	// Change the state of A and B.
	Particle* pA = system.storage->lookupLocalParticle(A);
	Particle* pB = system.storage->lookupLocalParticle(B);
	if (pA==NULL) {
	} else {
	  pA->setState(pA->getState()+deltaA);
	}
	if (pB==NULL) {
	} else {
	  pB->setState(pB->getState()+deltaB);
	}
	// Add a bond
	fpl->add(A, B);
      }

      LOG4ESPP_INFO(theLogger, "Leaving applyAR");

    }

    /****************************************************
     ** REGISTRATION WITH PYTHON
     ****************************************************/

    void AssociationReaction::registerPython() {
      using namespace espressopp::python;
      class_<AssociationReaction, shared_ptr<AssociationReaction>, bases<Extension> >
        ("integrator_AssociationReaction", init<shared_ptr<System>, shared_ptr<VerletList>, shared_ptr<FixedPairList>, shared_ptr<DomainDecomposition> >())
        .def("connect", &AssociationReaction::connect)
        .def("disconnect", &AssociationReaction::disconnect)
        .add_property("rate", &AssociationReaction::getRate, &AssociationReaction::setRate)
        .add_property("cutoff", &AssociationReaction::getCutoff, &AssociationReaction::setCutoff)
        .add_property("typeA", &AssociationReaction::getTypeA, &AssociationReaction::setTypeA)
        .add_property("typeB", &AssociationReaction::getTypeB, &AssociationReaction::setTypeB)
        .add_property("deltaA", &AssociationReaction::getDeltaA, &AssociationReaction::setDeltaA)
        .add_property("deltaB", &AssociationReaction::getDeltaB, &AssociationReaction::setDeltaB)
        .add_property("stateAMin", &AssociationReaction::getStateAMin, &AssociationReaction::setStateAMin)
        .add_property("interval", &AssociationReaction::getInterval, &AssociationReaction::setInterval)
        ;
    }
  }
}

