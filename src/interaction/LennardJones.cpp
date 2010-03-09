#include "LennardJones.hpp"

namespace espresso {
  namespace interaction {
    void LennardJonesFunction::preset() {
      real sig2 = sigma * sigma;
      real sig6 = sig2 * sig2 * sig2;
      ff1 = 48.0 * epsilon * sig6 * sig6;
      ff2 = 24.0 * epsilon * sig6;
      ef1 =  4.0 * epsilon * sig6 * sig6;
      ef2 =  4.0 * epsilon * sig6;
    }

    void LennardJonesFunction::setAutoShift() {
      real ratio = sigma / getCutoff();
      real rat2 = ratio * ratio;
      real rat6 = rat2 * rat2 * rat2;
      shift = 4.0 * epsilon * rat6 * (rat6 - 1.0);
    }

    void VerletListLennardJones::setAutoShift(int type1, int type2)
    {
      getFunction(type1, type2).setAutoShift();
    }
  }
}
