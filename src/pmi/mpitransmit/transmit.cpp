#include "pmi/types.hpp"
#include "pmi/transmit.hpp"
#include "pmi/exceptions.hpp"
#include "pmi/functions.hpp"
#include "pmi/worker_internal.hpp"

#include <vector>
#include <sstream>

#include <mpi.h>

using namespace std;

LOG4ESPP_LOGGER(mpilogger, "pmi.mpitransmit");

namespace pmi {
  namespace transmit {
    //////////////////////////////////////////////////
    // Status definition
    //////////////////////////////////////////////////
    const unsigned short STATUS_OK = 0;
    const unsigned short STATUS_USER_ERROR = 1;
    const unsigned short STATUS_INTERNAL_ERROR = 2;
    const unsigned short STATUS_OTHER_ERROR = 3;
    

    //////////////////////////////////////////////////
    // Command definition
    //////////////////////////////////////////////////
    const unsigned int CMD_END = 0;
    const unsigned int CMD_ASSOC_CLASS = 1;
    const unsigned int CMD_ASSOC_METHOD = 2;
    const unsigned int CMD_CREATE = 3;
    const unsigned int CMD_INVOKE = 4;
    const unsigned int CMD_DESTROY = 5;
  

    //////////////////////////////////////////////////
    // Basic function
    //////////////////////////////////////////////////
    WorkerIdType
    getWorkerId() {
      return MPI::COMM_WORLD.Get_rank();
    }

    //////////////////////////////////////////////////
    // Broadcasts
    //////////////////////////////////////////////////
    void
    broadcastCommand(const unsigned int command, 
		     const IdType classId, 
		     const IdType methodId, 
		     const IdType objectId) {
      unsigned int msg[4];
      msg[0] = command;
      msg[1] = classId;
      msg[2] = methodId;
      msg[3] = objectId;
    
      LOG4ESPP_INFO(mpilogger, "Controller broadcasts command ("	\
		    << msg[0] << ", "				\
		    << msg[1] << ", "				\
		    << msg[2] << ", "				\
		    << msg[3] << ")."				\
		    );
      MPI::COMM_WORLD.Bcast(msg, 4, MPI::UNSIGNED, getControllerMPIRank());
      LOG4ESPP_DEBUG(mpilogger, "Controller finished broadcast.");
    }
  
  
    void
    broadcastName(const string &name, const unsigned int length) {
      // broadcast the string
      LOG4ESPP_INFO(mpilogger, "Controller broadcasts name \"" << name << "\".");
      MPI::COMM_WORLD.Bcast(const_cast<char*>(name.c_str()), length, MPI::CHAR, getControllerMPIRank());
      LOG4ESPP_DEBUG(mpilogger, "Controller finished broadcast.");
    }
  
    void
    broadcastPMIObjectId(const IdType objectId) {
      unsigned int oid = objectId;
    
      LOG4ESPP_INFO(mpilogger, "Controller broadcasts object (OID = "	\
		    << oid << ")."					\
		    );
      MPI::COMM_WORLD.Bcast(&oid, 1, MPI::UNSIGNED, getControllerMPIRank());
      LOG4ESPP_DEBUG(mpilogger, "Controller finished broadcast.");
    }
    
#ifndef PMI_OPTIMIZE
    void 
    gatherStatus() {
      int size = MPI::COMM_WORLD.Get_size();
      unsigned short allStatus[size];
      LOG4ESPP_DEBUG(mpilogger, "Controller gathers status of workers.");
      MPI::COMM_WORLD.Gather(allStatus, 1, MPI::UNSIGNED_SHORT, 
			     allStatus, 1, MPI::UNSIGNED_SHORT, getControllerMPIRank());
      LOG4ESPP_DEBUG(mpilogger, "Controller gathered status of workers.");
      for (int i = 1; i < MPI::COMM_WORLD.Get_size(); i++) {
	if (allStatus[i] != STATUS_OK) {
	  LOG4ESPP_DEBUG(mpilogger, "Controller got status "	\
			 << allStatus[i] << " from worker "	\
			 << i << ".");
	  MPI::Status mpiStatus;
	  MPI::COMM_WORLD.Probe(i, 99, mpiStatus);
	  int length = mpiStatus.Get_count(MPI::CHAR);
	  LOG4ESPP_DEBUG(mpilogger,					\
			 "Controller waits to get a message of length "	\
			 << length << " from worker " << i << ".");
	  char s[length];
	  MPI::COMM_WORLD.Recv(s, length, MPI::CHAR, i, 99);
	  string what = s;
	  switch (allStatus[i]) {
	  case STATUS_USER_ERROR:
	    PMI_THROW_USER_ERROR(what);
	  case STATUS_INTERNAL_ERROR:
	    PMI_THROW_INTL_ERROR(what);
	  case STATUS_OTHER_ERROR:
	    PMI_THROW_INTL_ERROR("other error: " << what);
	  default:
	    PMI_THROW_INTL_ERROR("unknown status: " << what);
	  }
	}
      }
      LOG4ESPP_INFO(mpilogger, "Controller gathered all OK from workers.");
    }
#endif

    void 
    endWorkers() {
      broadcastCommand(CMD_END, NOT_DEFINED, NOT_DEFINED, NOT_DEFINED);
    }
    
    void 
    associateClass(const string &name, const IdType id) {
      unsigned int l = name.length()+1;
      broadcastCommand(CMD_ASSOC_CLASS, id, NOT_DEFINED, l);
      broadcastName(name, l);
    }
    

    void 
    associateMethod(const string &name, const IdType id) {
      int l = name.length()+1;
      broadcastCommand(CMD_ASSOC_METHOD, NOT_DEFINED, id, l);
      broadcastName(name, l);
    }
    
    
    void 
    create(const IdType classId, const IdType objectId) {
      broadcastCommand(CMD_CREATE, classId, NOT_DEFINED, objectId);
    }
    
    void 
    invoke(const IdType classId, const IdType methodId, 
	   const IdType objectId) {
      broadcastCommand(CMD_INVOKE, classId, methodId, objectId);
    }

    void 
    destroy(const IdType classId, const IdType objectId) {
      broadcastCommand(CMD_DESTROY, classId, NOT_DEFINED, objectId);
    }
    
    //////////////////////////////////////////////////
    // Receive
    //////////////////////////////////////////////////
    
    const unsigned int*
    receiveCommand() {
      static unsigned int msg[4];
      LOG4ESPP_DEBUG(mpilogger, "Worker " << getWorkerId()	\
		     << " waits to receive a command.");
      
      MPI::COMM_WORLD.Bcast(msg, 4, MPI::UNSIGNED, getControllerMPIRank());
      LOG4ESPP_INFO(mpilogger, "Worker " << getWorkerId()		\
		    << " received command ("				\
		    << msg[0] << ", "					\
		    << msg[1] << ", "					\
		    << msg[2] << ", "					\
		    << msg[3] << ")."					\
		    );
      return msg;
    }

    string
    receiveName(const unsigned int length) {
      char s[length];
      LOG4ESPP_DEBUG(mpilogger, printWorkerId()			\
		     << "waits to receive a name of length "	\
		     << length << ".");
      MPI::COMM_WORLD.Bcast(&s, length, MPI::CHAR, getControllerMPIRank());
      string name = s;
      LOG4ESPP_INFO(mpilogger, printWorkerId()			\
		    << "received name \"" << name << "\".");
      return name;
    }
    
    IdType
    receivePMIObjectId() {
      unsigned int oid;
      
      LOG4ESPP_DEBUG(mpilogger, "Worker " << getWorkerId()	\
		     << " waits to receive an object id.");
      
      MPI::COMM_WORLD.Bcast(&oid, 1, MPI::UNSIGNED, getControllerMPIRank());
      LOG4ESPP_INFO(mpilogger, "Worker " << getWorkerId()		\
		    << " received object id "				\
		    << oid << "."					\
		    );
      
      return oid;
    }
    
#ifndef PMI_OPTIMIZE
    void
    reportOk() {
      LOG4ESPP_INFO(mpilogger, printWorkerId()	\
		    << "reports status OK.");
      MPI::COMM_WORLD.Gather(&STATUS_OK, 1, MPI::UNSIGNED_SHORT, 
			     0, 0, MPI::UNSIGNED_SHORT, getControllerMPIRank());
    }
    
    void
    reportError(unsigned char status, const string &what) {
      LOG4ESPP_INFO(mpilogger, printWorkerId()				\
		    << "reports error status " << status << ".");
      MPI::COMM_WORLD.Gather(&status, 1, MPI::UNSIGNED_SHORT, 
			     0, 0, MPI::UNSIGNED_SHORT, getControllerMPIRank());
      LOG4ESPP_INFO(mpilogger, printWorkerId()			\
		    << "sends error message \"" << what << "\".");
      MPI::COMM_WORLD.Send(what.c_str(), what.length()+1, MPI::CHAR, getControllerMPIRank(), 99);
    }

    void reportInternalError(const string &what) {
      reportError(STATUS_INTERNAL_ERROR, what);
    }

    void reportUserError(const string &what) {
      reportError(STATUS_USER_ERROR, what);
    }
#endif
    
    bool
    handleNext() {
      const unsigned int *msg = receiveCommand();
      static string name;
      
      LOG4ESPP_DEBUG(mpilogger, "Worker " << getWorkerId()		\
		     << " handles command ("				\
		     << msg[0] << ", "					\
		     << msg[1] << ", "					\
		     << msg[2] << ", "					\
		     << msg[3] << ")"					\
		     );
      
      switch (msg[0]) {
      case CMD_END:
#ifndef PMI_OPTIMIZE
	reportOk();
#endif
	break;
      case CMD_ASSOC_CLASS:
	name = receiveName(msg[3]);
	worker::associateClass(name, msg[1]);
	break;
      case CMD_ASSOC_METHOD:
	name = receiveName(msg[3]);
	worker::associateMethod(name, msg[2]);
	break;
      case CMD_CREATE:
    worker::create(msg[1], msg[3]);
    break;
      case CMD_INVOKE:
    worker::invoke(msg[1], msg[2], msg[3]);
    break;
      case CMD_DESTROY:
	worker::destroy(msg[1], msg[3]);
	break;
      default:
	ostringstream ost; 
	ost << printWorkerId()
	    << " cannot handle unknown command code "
	    << msg[0] << ".";
	LOG4ESPP_FATAL(logger, ost.str());
#ifndef PMI_OPTIMIZE
	pmi::transmit::reportInternalError(ost.str());
#endif
      }
      
      return msg[0] != CMD_END;
    }
  }
}
