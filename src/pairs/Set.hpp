#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include "types.hpp"
#include "Computer.hpp"

namespace espresso {
  namespace storage {
    class Storage;
  }

  namespace pairs {
    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;
      
      Set() {}

      virtual ~Set() {}

      /** Apply computer to all pairs of this set. Call prepare() and
	  finalize().
       */
      virtual bool foreachPair(Computer& computer);
      virtual bool foreachPair(Computer::SelfPtr computer);

      /** approximate looping over all pairs of this set with a given
	  maximal distance. The pairs should always belong to this
	  pair set, but can exceed the maximal distance. A set should
	  try to avoid as many as possible pairs with a bigger
	  distance, however, for performance reasons.  This is not and
	  should not be exported to Python. */
      virtual bool enclForeachPairWithin(Computer &computer);

      /// storage that the lhs particles of the computer come from
      virtual shared_ptr< storage::Storage > getLeftStorage() = 0;
      /// storage that the rhs particles of the computer come from
      virtual shared_ptr< storage::Storage > getRightStorage() = 0;

      static void registerPython();

    protected:
      virtual bool foreachPairApply(Computer &computer) = 0;
    };
  }
}

#endif
