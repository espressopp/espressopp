#ifndef _PMI_HPP
#define _PMI_HPP

#include "acconfig.hpp"

#ifdef HAVE_MPI
#include <boost/mpi.hpp>
#include "pmi/pmi.hpp"

#define PMI_DECL_SETTER(method, type, arg)	\
  virtual void method(type arg);		\
  virtual void method##Worker();		\
  virtual void method##Local(type arg);

#define PMI_DEFINE_SETTER(_namespace, _class, method, type, arg)		\
  void _namespace::_class::method(type arg) {				\
    pmiObject.invoke<&_namespace::_class::method##Worker>();		\
    boost::mpi::communicator world;					\
    boost::mpi::broadcast(world, arg, pmi::getControllerMPIRank());	\
    method##Local(arg);							\
  }									\
									\
  void _namespace::_class::method##Worker() {				\
    type arg;								\
    boost::mpi::communicator world;					\
    boost::mpi::broadcast(world, arg, pmi::getControllerMPIRank());	\
    method##Local(arg);							\
  }									\
  PMI_REGISTER_METHOD(_namespace::_class, method##Worker); \
									\
  void _namespace::_class::method##Local(type arg)


#else // HAVE_MPI

#define PMI_REGISTER_CLASS(class)
#define PMI_REGISTER_METHOD(class, method)

#define PMI_DECL_SETTER(method, type, arg)	\
  virtual void method(type arg);
#define PMI_DEFINE_SETTER(_namespace, _class, method, type, arg)	\
  void _namespace::_class::method(type arg)


#endif // HAVE MPI

#endif
