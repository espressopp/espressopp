#include "espresso_common.hpp"

void initPythonEspresso()
{
  // the controller: register with python
  espresso::hello::registerPython();
}
