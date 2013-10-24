// ESPP_CLASS
#ifndef _ANALYSIS_STATICSTRUCTF_HPP
#define _ANALYSIS_STATICSTRUCTF_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "python.hpp"

namespace espresso {
    namespace analysis {

        /** Class to compute the static structure function of the system. */
        class StaticStructF : public Observable {
        public:
            StaticStructF(shared_ptr< System > system) : Observable(system) {
            }
            ~StaticStructF() {
            }
            virtual real compute() const;          
            virtual python::list computeArray(int nqx, int nqy, int nqz, real bin_factor) const;
            static void registerPython();
        };
    }
}

#endif
