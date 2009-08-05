#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include "types.hpp"
#include "Computer.hpp"
#include "storage/Storage.hpp"

namespace espresso {
  namespace pairs {
    class ForeachBreak: public std::exception {};

    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;
      
      Set(storage::Storage::SelfPtr _storage1,
	  storage::Storage::SelfPtr _storage2);

      virtual ~Set() {}

      /** Apply computer to all pairs of this set. Call prepare() and
	  finalize().
       */
      virtual bool foreach(Computer& computer);
      virtual bool foreach(Computer::SelfPtr computer);

      static void registerPython();

    protected:
      virtual bool foreachApply(Computer &computer) = 0;

      storage::Storage::SelfPtr storage1;
      storage::Storage::SelfPtr storage2;
    };
  }
}

#endif
