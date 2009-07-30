#ifndef _PARTICLES_SUBSET_HPP
#define _PARTICLES_SUBSET_HPP
#include "particles/Set.hpp"
#include "storage/Storage.hpp"

namespace espresso {
  namespace particles {
    class Subset : public Set {
    public:
      typedef shared_ptr< Subset > SelfPtr;
      
      Subset(const Set::SelfPtr superset);
      
      virtual storage::Storage::SelfPtr getStorage();

      static void registerPython();
    private:
      storage::Storage::SelfPtr storage;
    };
  }
}

#endif
