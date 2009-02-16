#ifndef _PMI_FUNCTIONS_HPP
#define _PMI_FUNCTIONS_HPP
/** @file Declares user functions of PMI. */

#include "pmi/types.hpp"
#include "pmi/exceptions.hpp"
#include "pmi/transmit.hpp"
#include "pmi/internal_func.hpp"

namespace pmi {
  /** Returns the MPI rank of the controller. Use this in broadcast
      and gather messages in MPI. */
  unsigned int getControllerMPIRank();

  /** Returns the worker id of the executing process. */
  WorkerIdType getWorkerId();

  /** Tests whether the executing process is the controller.  @return
      true, if the executing process is the controller, otherwise
      false. */
  bool isController();

  /** Tests whether the executing process is a worker. 
      @return true, if the executing process is a worker,
      otherwise false. */
  bool isWorker();

  /** Pretty prints the worker id. 
      @return a string that represents the worker id.
   */
  std::string printWorkerId();


  bool mainLoop();
  void endWorkers();

  /** Tests whether the workers are still active (i.e. whether
      endWorkers() has been called or not. 
      @return true, if the workers are still running. 
  */
  bool isWorkersActive();


  // Broadcast an object (call this only on the controller)
  template < class T >
  void 
  broadcastPMIObject(T& obj) {
    if (isController()) {
      pmi::transmit::broadcastPMIObjectId(obj.getObjectId());
    } else
      PMI_THROW_USER_ERROR("Called broadcastObject on worker!");
  }

  template < class T >
  T* getPMIObjectPtr(const IdType oid) {
    return static_cast<T*>(objects[oid]);
  }

  template < class T >
  T& receivePMIObject() {
    if (isController())
      PMI_THROW_USER_ERROR("Called receiveObject on controller!")
    else {
      IdType oid = pmi::transmit::receivePMIObjectId();
      return *getPMIObjectPtr<T*>(oid);
    }
  }

  template < class T >
  T* receivePMIObjectPtr() {
    if (isController())
      PMI_THROW_USER_ERROR("Called receiveObjectPtr on controller!")
    else {
      IdType oid = pmi::transmit::receivePMIObjectId();
      return getPMIObjectPtr<T>(oid);
    }
  }
  
  // Broadcast an object (SPMD style)
  // On the controller: broadcasts the object pointed to by ptr
  // On the worker: receive an object
  template < class T >
  void broadcastPMIObject(T*& ptr) {
    if (isWorker()) {
      ptr = receivePMIObjectPtr<T*>();
    } else {
      broadcastPMIObject<T>(*ptr);
    }
  }


}
#endif
