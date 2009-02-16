#include "pmi/worker_internal.hpp"

#include <vector>
#include "pmi/exceptions.hpp"
#include "pmi/transmit.hpp"
#include "pmi/functions.hpp"

using namespace pmi;
using namespace std;

// Define some simplified error reporting macros
#ifndef PMI_OPTIMIZE
#define PMI_REPORT_INTL_ERROR(output)		\
  {						\
    ostringstream ost; ost << output;		\
    LOG4ESPP_FATAL(logger, ost.str());					\
    pmi::transmit::reportInternalError(ost.str()); \
  }

#define PMI_REPORT_USER_ERROR(output)		\
  {						\
    ostringstream ost; ost << output;		\
    LOG4ESPP_FATAL(logger, ost.str());					\
    pmi::transmit::reportUserError(ost.str()); \
  }

#define PMI_REPORT_OK pmi::transmit::reportOk()

#else

#define PMI_REPORT_INTL_ERROR(output)		\
  {						\
    ostringstream ost; ost << output;		\
    LOG4ESPP_FATAL(logger, ost.str());					\
  }

#define PMI_REPORT_USER_ERROR(output)		\
  {						\
    ostringstream ost; ost << output;		\
    LOG4ESPP_FATAL(logger, ost.str());					\
  }

#define PMI_REPORT_OK

#endif

namespace pmi {
  namespace worker {
    vector<ConstructorCallerType> constructorCallers;
    vector<MethodCallerType> methodCallers;
    vector<DestructorCallerType> destructorCallers;

    void associateClass(const string &name, const IdType id) {
      if (constructorCallersByName().find(name) == 
	  constructorCallersByName().end())
	PMI_REPORT_USER_ERROR(printWorkerId()				\
			      << "has not registered class \"" << name	\
			      << "\" (constructor undefined).");
#ifndef PMI_OPTIMIZE
      if (destructorCallersByName().find(name) == destructorCallersByName().end())
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
			      << "has not registered class \"" << name	\
		       << "\" (constructur defined, but destructor undefined).");
      // check whether the vector has the right size
      if (constructorCallers.size() != id)
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "has " << constructorCallers.size()		\
		       << " associated constructors, but received id "	\
		       << id << " as the next id.");
      // check whether the vector has the right size
      if (destructorCallers.size() != id)
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "has " << destructorCallers.size()		\
		       << " associated destructors, but received id "	\
		       << id << " as the next id.");
      PMI_REPORT_OK;
#endif
      LOG4ESPP_INFO(logger, printWorkerId()			\
		    << "associates class \""			\
		    << name << "\" to class id " << id << ".");
      
      ConstructorCallerType cc = constructorCallersByName()[name];
      constructorCallers.push_back(cc);
      constructorCallersByName().erase(name);
      DestructorCallerType dc = destructorCallersByName()[name];
      destructorCallers.push_back(dc); 
      destructorCallersByName().erase(name);
    }
    
    void associateMethod(const string &name, const IdType id) {
      if (methodCallersByName().find(name) == methodCallersByName().end())
	PMI_REPORT_USER_ERROR(printWorkerId()				\
			      << "has not registered method \"" << name	\
			      << "\" (methodCaller undefined).");
#ifndef PMI_OPTIMIZE
      // check whether the vector has the right size
      if (methodCallers.size() != id)
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "has " << methodCallers.size()		\
		       << " associated methods, but received id "	\
		       << id << " as the next id.");
      PMI_REPORT_OK;
#endif
      LOG4ESPP_INFO(logger, printWorkerId()			\
		    << "associates method \""				\
		    << name << "\" to method id " << id << ".");
      
      MethodCallerType mc = methodCallersByName()[name];
      methodCallers.push_back(mc);
      methodCallersByName().erase(name);
    }
    
    void create(const IdType classId, 
		const IdType objectId) {
#ifndef PMI_OPTIMIZE
      // check whether the class has an associated constructorCaller
      if (classId >= constructorCallers.size())
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "does not have a constructorCaller for class id " \
		       << classId << ".");
      // check whether the object is already defined
      if (objectId < objects.size() && objects[objectId] != NULL) {
	PMI_REPORT_INTL_ERROR(printWorkerId()				\
		       << "has object id "			\
		       << objectId << " already defined.");
      } else if (objectId >= objects.size()+1) {
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << ": object id " << objectId			\
		       << " overshoots (size " << objects.size()	\
		       << ").");
      }
      PMI_REPORT_OK;
#endif
      LOG4ESPP_INFO(logger, printWorkerId()				\
		    << "creates an instance of class id "		\
		    << classId << ", object id is " << objectId << ".");
      
      // create object
      ConstructorCallerType cc = constructorCallers[classId];
      if (objectId == objects.size())
	objects.push_back((*cc)());
      else
	objects[objectId] = (*cc)();
    }
    
    void invoke(const IdType classId, 
		const IdType methodId, 
		const IdType objectId) {
#ifndef PMI_OPTIMIZE
      // check whether the method has an associated methodCaller
      if (methodId >= methodCallers.size())
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "does not have a methodCaller for method id " \
		       << methodId << ".");
      if (objectId >= objects.size() || objects[objectId] == NULL)
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "does not have an initialized object at object id " \
		       << objectId << " (invoke).");
      // TODO: check if object is of the right class
      PMI_REPORT_OK;
#endif
      // call method
      LOG4ESPP_INFO(logger, printWorkerId()				\
		    << "invokes method id " << methodId			\
		    << " of object id " << objectId			\
		    << " of class id " << classId << ".");
      
      void* voidPtr = objects[objectId];
      MethodCallerType mc = methodCallers[methodId];
      (*mc)(voidPtr);
    }
    
    void 
    destroy(const IdType classId,
	    const IdType objectId) {
#ifndef PMI_OPTIMIZE
      // check whether the class has an associated constructorCaller
      if (classId >= destructorCallers.size())
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "does not have a destructorCaller for class id " \
		       << classId << ".");
      // check if the object is defined
      if (objectId >= objects.size() || objects[objectId] == NULL)
	PMI_REPORT_INTL_ERROR(printWorkerId()					\
		       << "does not have an initialized object at object id " \
		       << objectId << " (destroy).");
      // TODO: check if object is of the right class
      PMI_REPORT_OK;
#endif
      // delete object
      LOG4ESPP_INFO(logger, printWorkerId()				\
		    << "destroys object id " << objectId		\
		    << " of class id " << classId << ".");
      
      void* voidPtr = objects[objectId];
      DestructorCallerType dc = destructorCallers[classId];
      (*dc)(voidPtr);
      objects[objectId] = NULL;
    }
  }
}
