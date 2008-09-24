#ifndef _PMI_PARALLELOBJECT_HPP
#define _PMI_PARALLELOBJECT_HPP
// A ParallelObject represents a parallel instance of a class.
// To be a able to use it, the parallel instance has to be first
// registered as a parallel class on the controller as well as on the
// workers.

#include "pmi/types.hpp"
#include "pmi/transmit.hpp"
#include "pmi/exceptions.hpp"
#include "pmi/controller.hpp"
#include "pmi/ParallelClass.hpp"
#include "pmi/ParallelMethod.hpp"

#ifdef CONTROLLER
using namespace pmi;

// macro to define proxy methods
#define PMI_PROXY_METHOD(method)					\
  void method() {							\
    ParallelObject<SubjectClass>::invoke<&SubjectClass::method>();	\
  }

namespace pmi {
  IdType generateObjectId();
  void freeObjectId(const IdType id);

  template < class T >
  class ParallelObject {
  public:
    typedef ParallelClass<T> PClass;
    typedef T SubjectClass;

  private:
    // the Id of the (parallel) instance
    IdType ID;
    // the instance running on the controller
    SubjectClass *objectPtr;

  public:
    ParallelObject() {
      if (isWorker())
	PMI_USER_ERROR(printWorkerId()					\
		       << "tries to create a parallel object of type \"" \
		       << PClass::getName() << "\".");
      if (!isWorkersActive())
	PMI_USER_ERROR("Controller tries to create a parallel object of type \"" \
		       << "\", but the workers have been terminated.");
      
      IdType classId = PClass::associate();

      // generate the object ID
      ID = generateObjectId();

      LOG4ESPP_INFO(logger, "Controller creates an instance of class \"" \
		    << PClass::getName() << "\" (class id " 		\
		    << PClass::getId() << "), object id is " << ID << ".");
      // broadcast invoke command
      transmit::create(classId, ID);

      // create the local instance
      objectPtr = new T();

#ifndef PMI_OPTIMIZE
      transmit::gatherStatus();
#endif
    }
      
    template < void (T::*method)() >
    void invoke() {
      if (isWorker())
	PMI_USER_ERROR(printWorkerId()					\
		       << "tries to invoke method \""			\
		       << (ParallelMethod<T,method>::getName())		\
		       << "\" of parallel object id " << ID		\
		       << " of class \"" << PClass::getName()		\
		       << "\".");

      if (!isWorkersActive())
	PMI_USER_ERROR("Controller tries to invoke method \""		\
		       << (ParallelMethod<T,method>::getName())		\
		       << "\" of parallel object id " << ID		\
		       << " of class \"" << PClass::getName()		\
		       << "\", but the workers have been terminated.");

      IdType methodId =
	ParallelMethod<T,method>::associate();

      IdType classId = PClass::getId();
      LOG4ESPP_INFO(logger, "Controller invokes method \""		\
		    << (ParallelMethod<T,method>::getName())		\
		    << "\" (method id " << methodId			\
		    << ") of object id " << ID				\
		    << " of class \"" << PClass::getName()		\
		    << "\" (class id " << classId << ").");
      // broadcast invoke command
      transmit::invoke(classId, methodId, ID);

      // invoke the method in the local objectPtr
      (objectPtr->*method)();

#ifndef PMI_OPTIMIZE
      transmit::gatherStatus();
#endif
    }

    ~ParallelObject() {
      if (isWorker())
	PMI_USER_ERROR(printWorkerId()					\
		       << "tries to destroy parallel object id " << ID	\
		       << " of class \"" << PClass::getName()		\
		       << "\".");

      if (isWorkersActive()) {
	IdType classId = PClass::getId();
	LOG4ESPP_INFO(logger, "Controller destroys object id " << ID	\
		      << " of class \"" << PClass::getName()		\
		      << "\" (class id " << classId << ").");
	transmit::destroy(classId, ID);
      } else {
	LOG4ESPP_DEBUG(logger, "Controller did not broadcast destroy message, as the workers are stopped.");
      }

      // destroy the local instance
      delete objectPtr;
      // free the objectId
      freeObjectId(ID);

#ifndef PMI_OPTIMIZE
      if (isWorkersActive()) {
	transmit::gatherStatus();
      } else {
	LOG4ESPP_DEBUG(logger, "Controller did not gather status after destroy, as the workers are stopped.");
      }
#endif
    }
  };
}

#endif /* CONTROLLER */
#endif /* _PMI_PARALLELOBJECT_HPP */
