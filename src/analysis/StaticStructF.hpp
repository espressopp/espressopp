// ESPP_CLASS
#ifndef _ANALYSIS_STATICSTRUCTF_HPP
#define _ANALYSIS_STATICSTRUCTF_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "python.hpp"

#include <sstream>//string convert for output files
#include <string>//should already be in sstream

#include <math.h>       //cos and ceil and sqrt
#include <algorithm>    // std::min
#include <functional>   // std::plus

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
            //computeArray usually is const, too!!!
            //virtual python::list computeArray(int nqx, int nqy, int nqz, real bin_factor);
            virtual python::list computeArray(int nqx, int nqy, int nqz, real bin_factor) const;
            static void registerPython();
            
        private:
//            void calculateWith(int index, ConfigurationPtr& c,
//                    double dq[], int nx, int ny, int nz, real& b_size,
//                    std::vector<int>& count_b, std::vector<real>& q_b, 
//                    std::vector<real>& s_b) const;           
            //void StaticStructF::calculateWith(int index);
 
            //static std::vector< std::vector< unsigned int > > makeOneIndex(int n);
            //needed for parallel test with right IDs (comment out other printVec):
            //real calcSFSum(int fnprocs, int fnum_part, 
            //                              Real3D fq, ConfigurationPtr fconfig) const;
            //void printVec(std::vector< std::vector<int> > v);      
            std::string convertInt(int number);
            void printVec(std::vector<real>& v) const;
        };
    }
}

#endif
