#include "espresso.hpp"

#include "hello.hpp"

void initPythonEspresso()
{
  // the controller: register with python
  hello::registerPython();
}
