#ifndef _PARTICLES_SUBSET_HPP
#define _PARTICLES_SUBSET_HPP
#include "particles/Set.hpp"
#include "particles/Storage.hpp"

namespace espresso {
  namespace particles {
    class Subset : public Set {
    public:
      typedef shared_ptr< Subset > SelfPtr;
      
      Subset(const Set::SelfPtr superset);
      
      virtual Storage::SelfPtr getStorage();

      static void registerPython();
    private:
      Storage::SelfPtr storage;
    };
  }
}

#endif
