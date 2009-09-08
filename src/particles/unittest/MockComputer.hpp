#ifndef _PARTICLES_UNITTEST_MOCKCOMPUTER_HPP
#define _PARTICLES_UNITTEST_MOCKCOMPUTER_HPP

#include "particles/Computer.hpp"
#include "storage/Storage.hpp"

#include <boost/unordered_set.hpp>

using namespace espresso;
using namespace espresso::storage;

struct MockComputer : particles::Computer {
  bool prepareCalled;
  bool finalizeCalled;
  int applyCalled;
  boost::unordered_set<ParticleHandle> parts;

  MockComputer() {
    prepareCalled = false;
    finalizeCalled = false;
    applyCalled = 0;
  }

  void prepare(Storage::SelfPtr) {
    prepareCalled = true;
  }

  void finalize() {
    finalizeCalled = true;
  }

  bool apply(ParticleHandle handle) {
    applyCalled++;
    parts.insert(handle);
    return true;
  }
};

#endif
