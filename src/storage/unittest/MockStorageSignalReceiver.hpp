#ifndef MOCKSTORAGERECEIVER_HPP
#define MOCKSTORAGERECEIVER_HPP

#include "storage/Storage.hpp"

using namespace espresso;
using namespace espresso::storage;

struct MockStorageSignalReceiver
  : public enable_shared_from_this< MockStorageSignalReceiver > {
  int handlesChangedCount;

  MockStorageSignalReceiver(): handlesChangedCount(0) {}
  void connect(Storage::SelfPtr store) {
    store->connections.add(
                           store->handlesChanged,
                           shared_from_this(),
                           &MockStorageSignalReceiver::receive);
  }

  void receive() { handlesChangedCount++; }
};

#endif
