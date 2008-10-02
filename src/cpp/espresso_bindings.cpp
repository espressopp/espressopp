#include "espresso.hpp"

void initPythonEspresso()
{
  // the controller: register with python
  hello::PHelloWorld::registerPython();
}
