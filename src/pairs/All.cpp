#include "All.hpp"
#include "esutil/virtual_functional.hpp"
#include "particles/Computer.hpp"
#include <python.hpp>

using namespace espresso;
using namespace espresso::bc;
using namespace espresso::pairs;
using namespace espresso::particles;

using namespace boost;

// some abbreviations

// Helper class 1

template <class ParticleComputer>
class Traverser1 
  : public ParticleComputer  {

  typedef typename ParticleComputer::argument_type Reference;
  typedef pairs::ComputerBase<Reference> PairComputer;
  
private:
  class Traverser2 
    : public ParticleComputer  {
    
  public:
    const BC& bc;
    
    const ConstPropertyHandle<ParticleId> id;
    const ConstPropertyHandle<Real3D> pos;

    const Reference pref1;
    const Real3D pos1;
    const ParticleId id1;

    PairComputer& pairComputer;

    Traverser2(const All* all,
	       const Reference pref,
	       PairComputer& _pairComputer
	       ) :
      bc(*(all->getBC().get())), 
      id(all->getSet()->getStorage()->getIdPropertyHandle()),
      pos(PropertyHandle<Real3D>(*(all->getCoordinateProperty()))),
      pref1(pref),
      pos1(pos[pref1]),
      id1(id[pref1]),
      pairComputer(_pairComputer)
    {  } 

    virtual void operator()(const Reference pref2) {
   
      if (id1 < id[pref2]) {
	Real3D dist = bc.getDist(pos1, pos[pref2]);
	pairComputer(dist, pref1, pref2);
      }
    }
  };

public:

  PairComputer& pairComputer;

  const All* all;

  Traverser1(const All* _all, PairComputer& _pairComputer) :
    pairComputer(_pairComputer),
    all(_all)
  {  }

  virtual void operator()(const Reference pref) {
    // printf ("Traverser1: will call traverser2\n");
    Traverser2 traverser2(all, pref, pairComputer);
    all->getSet()->foreach(traverser2);
  }
};

/*--------------------------------------------------------------------------
  All::~All()
  -------------------------------------------------------------------------- */

All::~All() {}

/*--------------------------------------------------------------------------
  All::All(boundary_conditions, particle_set)
  -------------------------------------------------------------------------- */

All::All(bc::PBC _bc,
         particles::PSet _set,
         PReal3DProperty _coordinates) :
  set(_set),
  bc(_bc),
  coordinates(_coordinates)
{ }

/*--------------------------------------------------------------------------
  All::foreach(PairComputer)
  -------------------------------------------------------------------------- */
void All::foreach(pairs::Computer& pairComputer) {
  // printf ("ParticlePairComputer non-const\n");
  Traverser1<particles::Computer> traverser1(this, pairComputer);;
  set->foreach(traverser1);
}
       
/*--------------------------------------------------------------------------
  All::foreach(ConstPairComputer)
  -------------------------------------------------------------------------- */
void All::foreach(pairs::ConstComputer& pairComputer) const {
  // printf ("ParticlePairComputer const\n");
  Traverser1<particles::ConstComputer> traverser1(this, pairComputer);;
  set->foreach(traverser1);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
All::registerPython() {
  using namespace boost;
  using namespace boost::python;

  class_<All, bases<Set> >
    ("pairs_All", 
     init< bc::PBC, particles::PSet, PReal3DProperty >())
    ;
}
