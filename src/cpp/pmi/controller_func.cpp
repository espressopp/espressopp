// Implements the broadcast
#include "pmi/controller_func.hpp"
#include "pmi/transmit.hpp"
#include "pmi/exceptions.hpp"

#include <vector>
#include <iostream>
using namespace std;
using namespace pmi;

bool workersActive = true;

bool 
pmi::isWorkersActive() { return workersActive; }

void 
pmi::endWorkers() {
  if (isWorker())
    PMI_USER_ERROR(printWorkerId() << "tries to end workers.");
  LOG4ESPP_INFO(logger, "Controller ends all workers.");
  transmit::endWorkers();
  workersActive = false;
#ifndef PMI_OPTIMIZE
  transmit::gatherStatus();
#endif
}
