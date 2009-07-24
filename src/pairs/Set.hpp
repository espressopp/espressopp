#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include "types.hpp"
#include "Computer.hpp"
#include "particles/Storage.hpp"

namespace espresso {
  namespace pairs {
    class ForeachBreak: public std::exception {};

    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;
      
      Set(particles::Storage::SelfPtr _storage1,
	  particles::Storage::SelfPtr _storage2);

      virtual ~Set() {}

      /** Apply computer to all pairs of this set. Call prepare() and
	  finalize().
       */
      virtual void foreach(Computer& computer);
      virtual void foreach(Computer::SelfPtr computer);

      virtual int bump() { return 1; }

      static void registerPython();

    protected:
      virtual void foreachApply(Computer &computer) = 0;

      particles::Storage::SelfPtr storage1;
      particles::Storage::SelfPtr storage2;
    };
  }
}

#endif
