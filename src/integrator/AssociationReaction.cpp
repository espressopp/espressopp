#include "python.hpp"
#include "AssociationReaction.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  namespace integrator {

    //using namespace espresso::iterator;

    AssociationReaction::AssociationReaction(
					     shared_ptr<System> system,
					     shared_ptr<VerletList> _verletList,
					     shared_ptr<FixedPairList> _fpl)
      :Extension(system), verletList(_verletList), fpl(_fpl) {

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
      LOG4ESPP_INFO(theLogger, "Perform AssociationReaction");
      //System& system = getSystemRef();
      //system.storage->updateGhostsV();

      // loop over VL pairs
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
	Particle &p1 = *it->first;
	Particle &p2 = *it->second;

	reactPair(p1, p2);
      }
    }

    /** For a given pair of particles, check if they meet the condition the
	reactive scheme. If it is so, the (A,B) pair is added to Alist.
    */
    void AssociationReaction::reactPair(Particle& p1, Particle& p2) {
      Real3D r = p1.position() - p2.position();
      real dist2 = r.sqr();
      bool found;
      found=false;
      if(dist2 < cutoff_sqr){
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
    void AssociationReaction::sendMultiMap(boost::unordered_multimap<int, int> &mm) {

      InBuffer inBuffer(getSystem()->comm);
      OutBuffer outBuffer(getSystem()->comm);

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

	real curCoordBoxL = getSystem()->bc->getBoxL()[coord];
    
	// lr loop: left right
	for (int lr = 0; lr < 2; ++lr) {
	  int dir         = 2 * coord + lr;
	  int oppositeDir = 2 * coord + (1 - lr);

	  if (nodeGrid.getGridSize(coord) == 1) {
	    LOG4ESPP_DEBUG(logger, "no communication");
	  }
	  else {
	    // prepare send and receive buffers
	    longint receiver, sender;
	    outBuffer.reset();
	    receiver = nodeGrid.getNodeNeighborIndex(dir);
	    sender = nodeGrid.getNodeNeighborIndex(oppositeDir);

	    // fill outBuffer from mm
	    outBuffer.write(mm.size());
	    for (reaclist::iterator it=mm.begin(); it!=mm.end(); it++) {
	      outBuffer.write(it->first);
	      outBuffer.write(it->second);
	    }

	    // exchange particles, odd-even rule
	    if (nodeGrid.getNodePosition(coord) % 2 == 0) {
	      outBuffer.send(receiver, AR_COMM_TAG);
	      inBuffer.recv(sender, AR_COMM_TAG);
	    } else {
	      inBuffer.recv(sender, AR_COMM_TAG);
	      outBuffer.send(receiver, AR_COMM_TAG);
	    }

	    // unpack received data
	    // add content of inBuffer to mm
	    int lengthA, Aidx, Bidx;
	    inBuffer.read(&lengthA);
	    for (int i=0; i<lengthA; i++) {
	      inBuffer.read(&Aidx);
	      inBuffer.read(&Bidx);
	      mm.insert(std::make_pair(Aidx, Bidx));
	    }
	  }
	}
      }
    }

    /** Collects, for each candidate B particle, all candidate A partners. Then
	picks a single A, deciding of the effective partners that are then stored in
	"partners".
    */
    void AssociationReaction::sortAndPickB() {

      // Collect reaction by B
      boost::unordered_multimap<int, int> Blist;
      boost::unordered_set<int> Bset;
      int A, B;

      // Create Blist, containing (B,A) pairs
      // Create Bset, containing each active B particle only once.
      for (reaclist::iterator it=Alist.begin(); it!=Alist.end(); it++) {
	A = it->first;
	B = it->second;
	Blist.insert(std::make_pair(B, A));
	Bset.append(B);
      }

      partners.erase();
      // For each active B, pick a partner
      for (boost::unordered_set<int>::iterator it=Bset.begin(); it!=Bset.end(); it++) {
	B = *it;
	int Bsize = Blist.count(B);
	if (Bsize>0) {
	  int pick = uni(Bsize);
	  std::pair<boost::unordered_multimap<int, int>::iterator,boost::unordered_multimap<int, int>::iterator> candidates;
	  candidates = Blist.equal_range(B);
	  int i=0;
	  for (boost::unordered_multimap<int, int>::iterator jt=candidates.first; jt!=candidates.second; jt++) {
	    if (i==pick) {
	      partners.add(jt->first, jt->second);
	      i++;
	      break
		}
	  }
	}
      }
    }

    /** Use the (A,B) list "partners" to add bonds and change the state of the
	particles accordingly.
    */
    void AssociationReaction::applyAR() {
      int A, B;
      for (boost::unordered_multimap<int, int>::iterator it=partners.begin(); it!=partners.end(); it++) {
	A = it->first;
	B = it->second;
	// Add a bond
	fpl->add(A, B);
	// Change the state of A and B.
	Particle &pA = lookupLocalParticle(A);
	Particle &pB = lookupLocalParticle(B);
	AAA->state() -= 1;
	BBB->state() += 1;
      }

    }

    /****************************************************
     ** REGISTRATION WITH PYTHON
     ****************************************************/

    void AssociationReaction::registerPython() {
      using namespace espresso::python;
      class_<AssociationReaction, shared_ptr<AssociationReaction>, bases<Extension> >
        ("integrator_AssociationReaction", init<shared_ptr<System>, shared_ptr<VerletList>, shared_ptr<FixedPairList> >())
        .def("connect", &AssociationReaction::connect)
        .def("disconnect", &AssociationReaction::disconnect)
        .add_property("rate", &AssociationReaction::getRate, &AssociationReaction::setRate)
        .add_property("cutoff", &AssociationReaction::getCutoff, &AssociationReaction::setCutoff)
        .add_property("typeA", &AssociationReaction::getTypeA, &AssociationReaction::setTypeA)
        .add_property("typeB", &AssociationReaction::getTypeB, &AssociationReaction::setTypeB)
        .add_property("deltaA", &AssociationReaction::getDeltaA, &AssociationReaction::setDeltaA)
        .add_property("deltaB", &AssociationReaction::getDeltaB, &AssociationReaction::setDeltaB)
        .add_property("stateAMin", &AssociationReaction::getStateAMin, &AssociationReaction::setStateAMin)
        ;
    }
  }
}

