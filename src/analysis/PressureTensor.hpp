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
      PressureTensor(shared_ptr< System > system) : Observable(system) {}
      ~PressureTensor() {}
      virtual real compute() const;
      virtual Tensor computeTensor() const;

      static void registerPython();
    };
  }
}

#endif
