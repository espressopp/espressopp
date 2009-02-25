#ifndef _PMI_PARALLELMETHOD_HPP
#define _PMI_PARALLELMETHOD_HPP
#include "pmi/ParallelClass.hpp"

// macro to register a method
#define PMI_REGISTER_METHOD(aClass, aMethod)				\
  template <>								\
  std::string pmi::ParallelMethod<aClass, &aClass::aMethod>::MNAME =	\
    pmi::ParallelMethod<aClass, &aClass::aMethod>			\
    ::registerMethod(#aClass "::" #aMethod "()");

namespace pmi { 
  IdType generateMethodId();

  template < class T, void (T::*method)() >
  class ParallelMethod {
  public:
    // store the name of the method
    static std::string MNAME;

    // store the Id of the method
    static IdType MID;

    // register the method
    // this is typically called statically
    static std::string registerMethod(const std::string &_name) {
      methodCallersByName()[_name] = methodCallerTemplate<T, method>;
      return _name;
    }

    static const std::string &getName() { return MNAME; }
    //    static IdType &getId() { return MID; }
    
    static IdType &associate() {
      if (MID == NOT_ASSOCIATED) {
	MID = generateMethodId();
	
	LOG4ESPP_INFO(logger, "Controller associates method \"" << MNAME << \
		      "\" to method id " << MID << ".");
	transmit::associateMethod(MNAME, MID);
#ifndef PMI_OPTIMIZE
	transmit::gatherStatus();
#endif
      }
      return MID;
    }
  };

  // Initialize MID
  template < class T, void (T::*method)() >
  IdType ParallelMethod<T, method>::MID = NOT_ASSOCIATED;

  // MNAME is not initialized: it has to be initialized when the method
  // is registered. Use PMI_REGISTER_METHOD for that purpose.

}

#endif
