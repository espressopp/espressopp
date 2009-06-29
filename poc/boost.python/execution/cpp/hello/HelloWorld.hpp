#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP

namespace hello {
  class HelloWorld {
  public:
    // print out the message
    void printMessage();
    // register the class in python
    static void registerPython();
    // functionality before python
    static void initBeforePython(int &argc, char **&argv);
  };

}

#endif
