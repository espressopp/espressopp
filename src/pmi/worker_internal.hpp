#ifndef _PMI_WORKER_INTERNAL_HPP
#define _PMI_WORKER_INTERNAL_HPP
#include "pmi/types.hpp"


namespace pmi {
  namespace worker {
    void associateClass(const std::string &name, const IdType id);
    void associateMethod(const std::string &name, const IdType id);
    
    void create(const IdType classId, const IdType objectId);
    void invoke(const IdType classId, 
		const IdType methodId, 
		const IdType objectId);
    void destroy(const IdType classId,
		 const IdType objectId);
  }
}

#endif /* _PMI_WORKER_INTERNAL_HPP */
