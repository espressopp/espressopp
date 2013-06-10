// ESPP_CLASS
#ifndef _ANALYSIS_ANALYSISBASE_HPP
#define _ANALYSIS_ANALYSISBASE_HPP
#include "types.hpp"
#include "python.hpp"
#include "SystemAccess.hpp"
#include <limits>

namespace espresso {
  namespace analysis {

    /** All quantities to be measured derive from this abstract base class and the
    * corresponding template.
    * The following functions have to be implemented by the derived class:
    * virtual ResultType computeRaw()
    * virtual python::list compute()
    * virtual python::list getAverageValue()
    * virtual void resetAverage()
    * virtual void updateAverage(ResultType res) = 0;
    * */
    class AnalysisBase : public SystemAccess {
    public:
      AnalysisBase(shared_ptr< System > system) : SystemAccess(system) {}
      virtual ~AnalysisBase() {}
      virtual void performMeasurement() = 0;
      virtual void reset() = 0;
      virtual python::list compute() = 0;
      virtual python::list getAverageValue() = 0;
      virtual int getNumberOfMeasurements() = 0;
      static void registerPython();
      static LOG4ESPP_DECL_LOGGER(logger);
    };

    template < class ResultType >
    class AnalysisBaseTemplate : public AnalysisBase {
    public:
      AnalysisBaseTemplate(shared_ptr< System > system) : AnalysisBase(system) {
    	  //this->reset();
      };
      virtual ~AnalysisBaseTemplate() {};
      virtual ResultType computeRaw() = 0;
      virtual void performMeasurement();
      virtual void reset();
      virtual python::list compute();
      virtual python::list getAverageValue();
      virtual void resetAverage();
      virtual void updateAverage(ResultType res);
      virtual int getNumberOfMeasurements();

      template <typename T>
        python::list computeTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == true ]);
      template <typename T>
        python::list computeTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == false ]);

      template <typename T>
        python::list getAverageValueTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == true ]);
      template <typename T>
        python::list getAverageValueTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == false ]);

      template <typename T>
        void resetAverageTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == true ]);
      template <typename T>
        void resetAverageTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == false ]);

      template <typename T>
        void updateAverageTypeSwitch(T res, char(*)[std::numeric_limits< T >::is_specialized == true]);
      template <typename T>
        void updateAverageTypeSwitch(T res, char(*)[std::numeric_limits< T >::is_specialized == false ]);

    protected:
      ResultType newAverage, lastAverage;
      ResultType newVariance, lastVariance;
      int nMeasurements;
    };

    /** This method will measure the observable and can be called from within the integrator
     * by using the ExtAnalyse extension. It can also be called from python.
     * It will also update the average value (typically also the variance) of the observable **/
    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::performMeasurement() {
    	ResultType res = computeRaw();
    	nMeasurements++;
    	updateAverage(res);
    }

    /** This method resets/initializes the measurement **/
    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::reset() {
    	nMeasurements = 0;
    	resetAverage();
    };

    /** return the number of measurements since the last reset **/
    template < class ResultType >
    inline int
    AnalysisBaseTemplate< ResultType >::getNumberOfMeasurements() {
      return nMeasurements;
    };

    /*
     * The following functions have to be overridden for non scalar ResultTypes
     * ( see PressureTensor.hpp as an example )
    */

    /*
     * python interface for computeRaw()
     */
    template < class ResultType >
    inline python::list
    AnalysisBaseTemplate< ResultType >::compute() {
      ResultType dummy;
      return computeTypeSwitch(dummy);
    }

    /*
     * python interface to get average and error bars
     */
    template < class ResultType >
    inline python::list
    AnalysisBaseTemplate< ResultType >::getAverageValue() {
      ResultType dummy;
      return getAverageValueTypeSwitch(dummy);
    }

    /*
     * reinitialize measurements
     */
    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::resetAverage() {
      ResultType dummy;
      resetAverageTypeSwitch(dummy);
      return;
    }

    /*
     * calculate new average and error bar after measurement
     */
    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::updateAverage(ResultType res) {
      updateAverageTypeSwitch(res);
      return;
    };

    /**
     *  Following the C++ SFINAE template pattern we make sure that the correct function
     *  is generated according to the condition:
     *
     *     std::numeric_limits< T >::is_specialized == true
     *
     *  (in general that means scalar integer and double types, I haven't testet it for complex.
     *  For all other types these functions have to be overridden -- see PressureTensor.hpp).
     *
    */

    template < class ResultType >
    template < class T >
    inline python::list
    AnalysisBaseTemplate< ResultType >::computeTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == true ] = 0) {
      python::list ret;
      ResultType res = computeRaw();
      ret.append(res);
      return ret;
    }

    template < class ResultType >
    template < class T >
    inline python::list
    AnalysisBaseTemplate< ResultType >::computeTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == false ] = 0) {
      python::list ret;
      return ret;
    }

    template < class ResultType >
    template < class T >
    inline python::list
    AnalysisBaseTemplate< ResultType >::getAverageValueTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == true ] = 0) {
      python::list ret;
      ResultType res = newAverage;
      ret.append(res);
      res = newVariance;
      ret.append(res);
      return ret;
    }

    template < class ResultType >
    template < class T >
    inline python::list
    AnalysisBaseTemplate< ResultType >::getAverageValueTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == false ] = 0) {
      python::list ret;
      return ret;
    }

    template < class ResultType >
    template < class T >
    inline void
    AnalysisBaseTemplate< ResultType >::resetAverageTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == true ] = 0) {
      newAverage   = 0;
      lastAverage  = 0;
      newVariance  = 0;
      lastVariance = 0;
      return;
    }

    template < class ResultType >
    template < class T >
    inline void
    AnalysisBaseTemplate< ResultType >::resetAverageTypeSwitch(T dummy, char(*)[std::numeric_limits< T >::is_specialized == false ] = 0) {
      return;
    }

    template < class ResultType >
    template < class T >
    inline void
    AnalysisBaseTemplate< ResultType >::updateAverageTypeSwitch(T res, char(*)[std::numeric_limits< T >::is_specialized == true ] = 0) {
    // compare e. g. Knuth TAOCP vol 2, 3rd edition, page 232
      if (nMeasurements == 1) {
          newAverage     = res;
          lastAverage    = newAverage;
      } else {
          newAverage  = lastAverage  + (res - lastAverage) / nMeasurements;
          newVariance = lastVariance + (res - lastAverage) * (res - newAverage);
          lastAverage = newAverage;
          lastVariance = newVariance;
      }
      return;
    };

    template < class ResultType >
    template < class T >
    inline void
    AnalysisBaseTemplate< ResultType >::updateAverageTypeSwitch(T res, char(*)[std::numeric_limits< T >::is_specialized == false] = 0) {
   	  return;
    };
  }
}
#endif
