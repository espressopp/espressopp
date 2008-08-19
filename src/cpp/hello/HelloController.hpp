#ifndef _HELLO_HELLOCONTROLLER_HPP
#define _HELLO_HELLOCONTROLLER_HPP

namespace hello {
  // class that calls a parallel HelloWorld class
  class HelloController {
  public:
    // print out the message
    void print();
    // register the controller class in python
    static void registerPython();
  };

}

#endif
