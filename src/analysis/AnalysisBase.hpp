// ESPP_CLASS
#ifndef _ANALYSIS_ANALYSISBASE_HPP
#define _ANALYSIS_ANALYSISBASE_HPP
#include "types.hpp"
#include "python.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace analysis {
    /** All quantities to be measured derive from this abstract base class. */

    class AnalysisBase : public SystemAccess {
    public:
      AnalysisBase(shared_ptr< System > system) : SystemAccess(system) {};
      virtual ~AnalysisBase() {};
      virtual void analyze() = 0;
      virtual python::list compute() = 0;
      virtual void reset() = 0;
      virtual int getNumberOfMeasurements() = 0;
      static void registerPython();
      static LOG4ESPP_DECL_LOGGER(logger);
    };

    template < class ResultType >
    class AnalysisBaseTemplate : public AnalysisBase {
    public:
      AnalysisBaseTemplate(shared_ptr< System > system) : AnalysisBase(system) {};
      virtual ~AnalysisBaseTemplate() {};
    public:
      virtual ResultType computeRaw() = 0;
      ResultType getAverage();
      ResultType getVariance();
      virtual void analyze();
      virtual python::list compute();
      virtual void reset();
      virtual int getNumberOfMeasurements();
    protected:
      int newAverage, lastAverage;
      int newVariance, lastVariance;
      int nMeasurements;
    };

    template < class ResultType >
    inline python::list
    AnalysisBaseTemplate< ResultType >::compute() {
        python::list ret;
        ResultType res = computeRaw();
        // BOOST_FOREACH(real value, result_real_vector) ret.append(value);
        ret.append(res);
        return ret;
    }


    template < class ResultType >
    inline ResultType
    AnalysisBaseTemplate< ResultType >::getAverage() {};

    template < class ResultType >
    inline ResultType
    AnalysisBaseTemplate< ResultType >::getVariance() {};

    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::reset() {};

    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::analyze() {
    	ResultType res = computeRaw();
    	return;
    };

    template < class ResultType >
    inline int
    AnalysisBaseTemplate< ResultType >::getNumberOfMeasurements() { return nMeasurements; };

  }
}

#endif
