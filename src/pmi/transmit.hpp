#ifndef _PMI_TRANSMIT_HPP
#define _PMI_TRANSMIT_HPP

#include "pmi/types.hpp"

namespace pmi {
  namespace transmit {
    // get the workerId of the executing process
    // corresponds to the MPI rank
    // 0 is the controller
    WorkerIdType getWorkerId();

    // The following functions should transmit the corresponding
    // commands to the workers 
    void endWorkers();

    void associateClass(const std::string &name, const IdType id);
    void associateMethod(const std::string &name, const IdType id);
    
    void create(const IdType classId,
		const IdType objectId);
    void invoke(const IdType classId, 
		const IdType methodId, 
		const IdType objectId);
    void destroy(const IdType classId,
		 const IdType objectId);
    void broadcastPMIObjectId(const IdType objectId);

#ifndef PMI_OPTIMIZE
    // collect the results of the last operation from all workers
    // check for failure
    void gatherStatus();

    void reportOk();
    void reportUserError(const std::string &what);
    void reportInternalError(const std::string &what);
#endif

    IdType receivePMIObjectId();

    // receive the next command from the controller and execute it
    // returns false if the stop worker command was received, true otherwise
    bool handleNext();
  }
}

#endif
