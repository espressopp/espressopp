#ifndef _PMI_PARALLELCLASS_HPP
#define _PMI_PARALLELCLASS_HPP
#include "pmi/types.hpp"
#include "pmi/transmit.hpp"
#include "pmi/exceptions.hpp"

#ifdef WORKER
#include "pmi/worker_func.hpp"
#endif

#include <iostream>

using namespace std;
using namespace pmi;

// macro to register a class
#define PMI_REGISTER_CLASS(aClass, name)					\
  template <>								\
  string pmi::ParallelClass<aClass>::NAME =				\
    pmi::ParallelClass<aClass>::registerClass(name);
  
namespace pmi {
  IdType generateClassId();

  // This class registers a parallel class with the Controller
  template < class T >
  class ParallelClass {
#ifdef CONTROLLER
    // store the name of the class
    static string NAME;

    // store the Id of the class
    static IdType ID;
#endif

  public:
    // register the class
    static const string &registerClass(const string &name) {
#ifdef WORKER
      // register constructorCaller with name
      constructorCallersByName()[name] = constructorCallerTemplate<T>;
      // register destructorCaller with name
      destructorCallersByName()[name] = destructorCallerTemplate<T>;
#endif
      return name;
    }

#ifdef CONTROLLER
    static const string &getName() { return NAME; }
    static IdType &getId() { return ID; }

    static IdType &associate() {
      if (ID == NOT_ASSOCIATED) {
	// associate class with ID
	ID = generateClassId();
	
	LOG4ESPP_INFO(logger, "Controller associates class \"" << NAME << \
		      "\" with class id " << ID << ".");
	transmit::associateClass(NAME, ID);
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
  template < class T >
  IdType ParallelClass<T>::ID = NOT_ASSOCIATED;
  // NAME is not initialized: it has to be initialized when the class
  // is registered. Use PMI_REGISTER_CLASS for that purpose.
#endif
}

#endif
