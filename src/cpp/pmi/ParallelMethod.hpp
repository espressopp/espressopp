#ifndef _PMI_PARALLELMETHOD_HPP
#define _PMI_PARALLELMETHOD_HPP
#include "pmi/types.hpp"
#include "pmi/transmit.hpp"
#include "pmi/exceptions.hpp"
#include "pmi/ParallelClass.hpp"

#ifdef WORKER
#include "pmi/worker_func.hpp"
#endif

using namespace std;
using namespace pmi;

// macro to register a class
#define PMI_REGISTER_METHOD(aClass, aMethod, name)			\
  template <>								\
  string pmi::ParallelMethod<aClass, &aClass::aMethod>::NAME =		\
    pmi::ParallelMethod<aClass, &aClass::aMethod>::registerMethod(name);

namespace pmi { 
  IdType generateMethodId();

  template < class T, void (T::*method)() >
  class ParallelMethod {
  public:
#ifdef CONTROLLER
    // store the name of the method
    static string NAME;

    // store the Id of the method
    static IdType ID;
#endif

    // register the method
    // this is typically called statically
    static string registerMethod(const string &_name) {
      string name = ParallelClass<T>::getName() + "::" +_name + "()";
#ifdef WORKER
      methodCallersByName()[name] = methodCallerTemplate<T, method>;
#endif
      return name;
    }

#ifdef CONTROLLER
    static const string &getName() { return NAME; }
    static IdType &getId() { return ID; }
    
    static IdType &associate() {
      if (ID == NOT_ASSOCIATED) {
	if (NAME == NOT_REGISTERED) {
	  LOG4ESPP_FATAL(logger, "Controller tried to associate a method that was not registered!");
	  throw UserError("Controller tried to associate a method that was not registered!");
	}
	ID = generateMethodId();
	
	LOG4ESPP_INFO(logger, "Controller associates method \"" << NAME << \
		      "\" to method id " << ID << ".");
	transmit::associateMethod(NAME, ID);
#ifndef PMI_OPTIMIZE
	transmit::gatherStatus();
#endif
      }
      return ID;
    }
#endif
  };

#ifdef CONTROLLER
  // Initialize ID
  template < class T, void (T::*method)() >
  IdType ParallelMethod<T, method>::ID = NOT_ASSOCIATED;
}
#endif

#endif
