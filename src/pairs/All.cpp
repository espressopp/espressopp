#include "All.hpp"

#include <exception>

#include "particles/Computer.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::bc;
using namespace espresso::pairs;
//using namespace espresso::particles;

using namespace boost;

All::All(bc::BC::SelfPtr _bc,
         particles::Set::SelfPtr _set,
         Property< Real3D >::SelfPtr _posProperty ) 
  : Set(_set->getStorage(), _set->getStorage()),
    set(_set), bc(_bc), posProperty(_posProperty)
{ }

namespace {
  class SameId {};

  template < class ApplyFunction, class ParticleHandle >
  class Traverser {
    typedef Traverser< ApplyFunction, ParticleHandle > Self;

    // the function to be applied to the pair
    const ApplyFunction &applyFunction;
    // the boundary conditions required to compute the posititon
    const bc::BC &bc;
    // the set to loop over
    particles::Set &set;

    particles::PropertyHandle< Real3D > pos;
    // the current particle in the first set
    ParticleHandle p1;

  public:
    Traverser(const ApplyFunction &_applyFunction,
	      const bc::BC::SelfPtr bcptr,
	      const particles::Set::SelfPtr setptr,
	      const Property< Real3D >::SelfPtr posProperty)
      : applyFunction(_applyFunction), 
	bc(*bcptr), 
	set(*setptr), 
	pos(*posProperty)
    {}

    void foreach() {
      particles::ApplyFunction af = 
	boost::bind(&Self::traverse1, this, _1);
      set.foreach(af);
    }

  private:
    void traverse1(const ParticleHandle _p1) {
      p1 = _p1;
      try {
	set.foreach(boost::bind(&Self::traverse2, this, _1));
      } catch (SameId) {}
    }
    
    void traverse2(const ParticleHandle p2) {
      // if a pair of a particle with itself turn up in the inner loop,
      // interrupt the inner loop and continue with the next element from
      // the outer loop
      if (p1 == p2) throw SameId();
      Real3D pos1 = pos[p1];
      Real3D pos2 = pos[p2];
      Real3D dist = bc.getDist(pos1, pos2);
      applyFunction(dist, p1, p2);
    }


  };

}

void All::foreach(ApplyFunction applyFunction) {
  Traverser< ApplyFunction, particles::ParticleHandle > 
    traverser(applyFunction, bc, set, posProperty);
  traverser.foreach();
}

void All::foreach(ConstApplyFunction applyFunction) const {
  Traverser< ConstApplyFunction, particles::ConstParticleHandle > 
    traverser(applyFunction, bc, set, posProperty);
  traverser.foreach();
}  

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
All::registerPython() {
  using namespace espresso::python;
  class_< All, bases< Set > >
    ("pairs_All", 
     init < bc::BC::SelfPtr, 
     particles::Set::SelfPtr, 
     Property< Real3D >::SelfPtr >())
    ;
}
