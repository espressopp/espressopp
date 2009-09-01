#ifndef SINGLE_NODE_HPP
#define SINGLE_NODE_HPP

#include "Storage.hpp"

namespace espresso {
  namespace storage {
    class SingleNode: public Storage {
    public:
      SingleNode(bc::BC::SelfPtr);
      virtual ParticleHandle getParticleHandle(ParticleId id) ;
      virtual ParticleHandle addParticle(ParticleId);
      virtual void deleteParticle(ParticleId);

      /// make this class available at Python
      static void registerPython();

    protected:
      virtual bool foreachApply(particles::Computer &);
      virtual esutil::TupleVector &getTupleVector();

    private:
      /// here the particle data is stored
      esutil::TupleVector particles;
    };
  }
}

#endif
