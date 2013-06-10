// ESPP_CLASS
#ifndef _ANALYSIS_OBSERVABLE_HPP
#define _ANALYSIS_OBSERVABLE_HPP
#include "python.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"
#include <vector>

namespace espresso {
  namespace analysis {
    /** All quantities to be measured derive from this abstract base class. */
    class Observable : public SystemAccess {
    public:
      /** Observable can be int, real, scalar or vector */
      enum result_types {old_format=-1, none=0, real_scalar=1, int_scalar=2, real_vector=3, int_vector=4};
      Observable(shared_ptr< System > system) : SystemAccess(system) {result_type=old_format; };
      virtual ~Observable() {};

    public:
      /** for compatibilty with old compute function only, will be deleted soon */
      virtual real compute() const { return 0.0; };
      /** returns observable of type real, used for Python and on C++ level*/
      virtual real compute_real() const { return 0.0;};
      /** returns observable of type int, used for Python and on C++ level*/
      virtual int compute_int() const { return 0; };
      /** computes vector of real values (e.g. pressure tensor, ...), used on C++ level only */
      virtual void compute_real_vector(){ return; };
      /** computes vector of integer values, used on C++ level only */
      virtual void compute_int_vector(){ return; };

      /** returns python list of real values (e.g. pressure tensor, ...), used on Python level*/
      virtual python::list compute_real_vector_python();
      /** returns python list of integer values, used on Python level*/
      virtual python::list compute_int_vector_python();

      /** returns the result type of the observable */
      // TODO at the moment it returns int instead, because it was causing an error
      // trying to convert result_types to python
      //result_types getResultType() { return result_type; };
      int getResultType() { return result_type; };

      static void registerPython();

     protected:
      result_types result_type;
      std::vector< real > result_real_vector;
      std::vector< int > result_int_vector;

      static LOG4ESPP_DECL_LOGGER(logger);

    };
  }
}

#endif
