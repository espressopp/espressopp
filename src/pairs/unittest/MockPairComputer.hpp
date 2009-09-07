#ifndef _PAIRS_UNITTEST_MOCKPAIRCOMPUTER_HPP
#define _PAIRS_UNITTEST_MOCKPAIRCOMPUTER_HPP

#include "pairs/Computer.hpp"
#include "storage/Storage.hpp"

#include <boost/unordered_set.hpp>

using namespace espresso;
using namespace espresso::storage;

struct MockPairComputer : pairs::Computer {
  Storage::SelfPtr storage;
  bool prepareCalled;
  bool finalizeCalled;
  int applyCalled;
  boost::unordered_set< std::pair<ParticleHandle, ParticleHandle> > pairs;

  MockPairComputer() {
    prepareCalled = false;
    finalizeCalled = false;
    applyCalled = 0;
  }

  void prepare(Storage::SelfPtr store, Storage::SelfPtr) {
    storage = store;
    prepareCalled = true;
  }

  void finalize() {
    finalizeCalled = true;
  }

  bool apply(const Real3D &,
             ParticleHandle handle1,
             ParticleHandle handle2
             ) {
    applyCalled++;
    pairs.insert(std::make_pair(handle1, handle2));
    return true;
  }
};

#endif
