#ifndef _MPIPO_WORKER_FUNC_HPP
#define _MPIPO_WORKER_FUNC_HPP

#include "pmi/types.hpp"

#ifdef WORKER
#include <map>

using namespace std;

namespace pmi {
  typedef void* (*ConstructorCallerType)();
  typedef void (*MethodCallerType)(void*);
  typedef void (*DestructorCallerType)(void*);

  // construct-on-first-use idioms
  map<string, ConstructorCallerType> &constructorCallersByName();
  map<string, MethodCallerType> &methodCallersByName();  
  map<string, DestructorCallerType> &destructorCallersByName();

  // templates for the different callers 
  template <class T>
  void*
  constructorCallerTemplate() {
    return static_cast<void*>(new T());
  }

  template <class T>
  void
  destructorCallerTemplate(void* voidPtr) {
    T* ptr = static_cast<T*>(voidPtr);
    delete ptr;
  }

  template <class T, class returnType, returnType (T::*methodPtr)()>
  void
  methodCallerTemplate(void *voidPtr) {
    T* ptr = static_cast<T*>(voidPtr);
    (ptr->*methodPtr)();
  }

  bool mainLoop();
}
#endif

#endif
