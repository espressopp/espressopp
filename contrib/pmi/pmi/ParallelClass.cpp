#include "ParallelClass.hpp"

namespace pmi {
  IdType generateClassId() {
    static IdType nextClassId = 0;
    nextClassId++;
    return nextClassId - 1;
  }
}
