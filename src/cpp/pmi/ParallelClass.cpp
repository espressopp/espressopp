#include "ParallelClass.hpp"
#include <set>

using namespace std;
using namespace pmi;

namespace pmi {
  IdType _generateClassId() {
    static IdType nextClassId = 0;
    nextClassId++;
    return nextClassId - 1;
  }

  // global object that manages the object IDs
  set<IdType> freeObjectIds;

  IdType 
  _generateObjectId() {
    static IdType nextObjectId = 0;
    if (freeObjectIds.empty()) {
      nextObjectId++;
      return nextObjectId - 1;
    } else {
      IdType id = *(freeObjectIds.begin());
      freeObjectIds.erase(id);
      return id;
    }
  }
  
  void 
  _freeObjectId(const IdType id) {
#ifndef PMI_OPTIMIZE
    if (freeObjectIds.find(id) != freeObjectIds.end())
      PMI_THROW_INTL_ERROR("Controller tried to free object id "	\
			   << id << " that is already free!");
#endif
    freeObjectIds.insert(id);
  }
}
