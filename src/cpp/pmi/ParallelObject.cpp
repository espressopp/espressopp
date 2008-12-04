#include "pmi/ParallelObject.hpp"

#ifdef CONTROLLER

#include <set>

using namespace std;
using namespace pmi;

set<IdType> freeObjectIds;

pmi::IdType 
pmi::generateObjectId() {
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
pmi::freeObjectId(const IdType id) {
#ifndef PMI_OPTIMIZE
  if (freeObjectIds.find(id) != freeObjectIds.end())
    PMI_INTL_ERROR("Controller tried to free object id "	\
		   << id << " that is already free!");
#endif
  freeObjectIds.insert(id);
}
#endif
