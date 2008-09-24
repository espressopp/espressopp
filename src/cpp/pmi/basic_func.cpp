#include "pmi/basic_func.hpp"
#include "pmi/transmit.hpp"
#include <sstream>

using namespace std;
using namespace pmi;

bool 
pmi::isController() {
  return getWorkerId() == CONTROLLER_ID;
}

bool 
pmi::isWorker() {
  return !isController();
}

// pretty print the worker id
string
pmi::printWorkerId() {
  if (isController())
    return "Controller ";
  else {
    ostringstream ost;
    ost << "Worker " << getWorkerId() << " ";
    return ost.str();
  }
}


