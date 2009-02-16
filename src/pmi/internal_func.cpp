#include "internal_func.hpp"

using namespace std;

namespace pmi {
  // These store the callers between registration and association
  // construct-on-first-use idiom
  map<string, ConstructorCallerType> &constructorCallersByName() {
    static map<string, ConstructorCallerType> callers;
    return callers;
  }
  
  // construct-on-first-use idiom
  map<string, MethodCallerType> &methodCallersByName() {
    static map<string, MethodCallerType> callers;
    return callers;
  }
  
  // construct-on-first-use idiom
  map<string, DestructorCallerType> &destructorCallersByName() {
    static map<string, DestructorCallerType> callers;
    return callers;
  }
  
  vector<void*> objects;
}  
