// ESPP_CLASS
#ifndef _ANALYSIS_PRESSURE_TENSOR_HPP
#define _ANALYSIS_PRESSURE_TENSOR_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "Tensor.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the pressure tensor. */
    class PressureTensor : public Observable {
    public:
      PressureTensor(shared_ptr< System > system) : Observable(system) {
        result_type=real_vector;
      }
      ~PressureTensor() {}
      
      virtual Tensor computeTensor() const;
      virtual void compute_real_vector();
      
      virtual python::list computeTensorIKz1(int, real) const;
      virtual Tensor computeTensorIKz2(real,real) const;

      static void registerPython();
    };
  }
}

#endif
