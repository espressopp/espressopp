#ifndef _PMI_PARALLELCLASS_HPP
#define _PMI_PARALLELCLASS_HPP
#include "pmi/functions.hpp"

#include <iostream>

// macro to register a class 
#define PMI_REGISTER_CLASS(aClass)				\
  template <>							\
  std::string pmi::ParallelClass<aClass>::CNAME =		\
    pmi::_registerClass<aClass>(#aClass);

// macro to create an spmd method
#define PMI_CREATE_SPMD_METHOD(name, class, method, object)	\
  void name() {							\
    object.invoke<&class::method>();				\
    this->method();						\
  }

namespace pmi {
  template < class T, void (T::*method)() >
  class ParallelMethod;

  // Functions used internally in the template
  // Don't use these functions!
  IdType _generateClassId();
  IdType _generateObjectId();
  void _freeObjectId(const IdType id);
  template <typename T>
  const std::string &_registerClass(const std::string &name) {
    // register constructorCaller with name
    constructorCallersByName()[name] = constructorCallerTemplate<T>;
    // register destructorCaller with name
    destructorCallersByName()[name] = destructorCallerTemplate<T>;
    return name;
  }

  // The ParallelClass template represents a parallel version of the
  // parameter class T (called the "subject class").
  // Each member function of the subject class can be invoked in
  // parallel via "invoke"
  template < class T >
  class ParallelClass {
  public:
    typedef T SubjectClass;

  private:
    // store the name of the class
    static std::string CNAME;

    // store the Id of the class
    static IdType CID;

    // store the Id of the object instance
    IdType OID;

  public:
    // returns the name of the class
    static const std::string &getName() { return CNAME; }
    IdType getObjectId() const { return OID; }

    // The constructor of the class: create an instance of the parallel object
    ParallelClass() {
      if (!isWorker()) {
	if (!isWorkersActive())
	  PMI_THROW_USER_ERROR("Controller tries to create a parallel object of type \"" \
			       << "\", but the workers have been terminated.");
	
	// Associate the class Id if not already done
	if (CID == NOT_ASSOCIATED) {
	  // associate class with ID
	  CID = _generateClassId();
	
	  LOG4ESPP_INFO(pmi::logger, "Controller associates class \"" << CNAME << \
			"\" with class id " << CID << ".");
	  transmit::associateClass(CNAME, CID);
#ifndef PMI_OPTIMIZE
	  transmit::gatherStatus();
#endif
	}
	
	// Generate the object ID
	OID = _generateObjectId();
	
	LOG4ESPP_INFO(logger, "Controller creates an instance of class \"" \
		      << CNAME << "\" (class id "			\
		      << CID << "), object id is " << OID << ".");
	// Broadcast creation command
      transmit::create(CID, OID);
      
#ifndef PMI_OPTIMIZE
      transmit::gatherStatus();
#endif
      }
    }

    template < void (SubjectClass::*method)() >
    void invoke() {
      const std::string &MNAME 
	= ParallelMethod<SubjectClass, method>::getName();
      
      if (isWorker())
	PMI_THROW_USER_ERROR(printWorkerId()				\
			      << "tries to invoke method \""		\
			      << MNAME					\
			      << "\" of parallel object id " << OID	\
			      << " of class \"" << CNAME		\
			      << "\".");

      if (!isWorkersActive())
	PMI_THROW_USER_ERROR("Controller tries to invoke method \""	\
			     << MNAME					\
			     << "\" of parallel object id " << OID	\
			     << " of class \"" << CNAME			\
			     << "\", but the workers have been terminated.");

      IdType MID =
	ParallelMethod<SubjectClass, method>::associate();

      LOG4ESPP_INFO(logger, "Controller invokes method \""		\
		    << MNAME					\
		    << "\" (method id " << MID				\
		    << ") of object id " << OID				\
		    << " of class \"" << CNAME				\
		    << "\" (class id " << CID << ").");
      // broadcast invoke command
      transmit::invoke(CID, MID, OID);
      
#ifndef PMI_OPTIMIZE
      transmit::gatherStatus();
#endif
    }

    ~ParallelClass() {
      if (!isWorker()) {
	if (isWorkersActive()) {
	  LOG4ESPP_INFO(logger, "Controller destroys object id " << OID	\
			<< " of class \"" << CNAME			\
			<< "\" (class id " << CID << ").");
	  transmit::destroy(CID, OID);
	} else {
	  LOG4ESPP_DEBUG(logger, "Controller did not broadcast destroy message, as the workers are stopped.");
	}
	
	// free the objectId
	_freeObjectId(OID);

#ifndef PMI_OPTIMIZE
	if (isWorkersActive()) {
	  transmit::gatherStatus();
	} else {
	  LOG4ESPP_DEBUG(logger, "Controller did not gather status after destroy, as the workers are stopped.");
	}
#endif
      }
    }
  };

  // Initialize class ID
  template < class T >
  IdType ParallelClass<T>::CID = NOT_ASSOCIATED;
  // CNAME is not initialized: it has to be initialized when the class
  // is registered. Use PMI_REGISTER_CLASS for that purpose.
}

#endif
