#include "ParallelMethod.hpp"

namespace pmi {
  IdType generateMethodId() {
    static IdType nextMethodId = 0;
    nextMethodId++;
    return nextMethodId - 1;
  }
}
