#ifndef _REALTENSOR_HPP
#define _REALTENSOR_HPP

#include "estypes.hpp"
#include "esutil/Tensor.hpp"

namespace espresso {
    /// make RealTensor globally accessible for convenience
    typedef esutil::Tensor<real> RealTensor;

    void registerPythonRealTensor();
}

#endif
