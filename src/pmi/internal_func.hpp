#ifndef _PMI_INTERNAL_FUNC_HPP
#define _PMI_INTERNAL_FUNC_HPP
/** @file Declares functions that should be used only internally by
    PMI. */

#include "pmi/types.hpp"
#include <map>
#include <vector>

namespace pmi {
  typedef void* (*ConstructorCallerType)();
  typedef void (*MethodCallerType)(void*);
  typedef void (*DestructorCallerType)(void*);

  // construct-on-first-use idioms
  std::map<std::string, ConstructorCallerType> &constructorCallersByName();
  std::map<std::string, MethodCallerType> &methodCallersByName();  
  std::map<std::string, DestructorCallerType> &destructorCallersByName();

  /** @internal
      Stores the object pointers of each object id. */
  extern std::vector<void*> objects;

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

  template <class T, void (T::*methodPtr)()>
  void
  methodCallerTemplate(void *voidPtr) {
    T* ptr = static_cast<T*>(voidPtr);
    (ptr->*methodPtr)();
  }
}

#endif
