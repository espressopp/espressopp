// ESPP_CLASS
#ifndef _ANALYSIS_PARTICLERADIUSDISTRIBUTION_HPP
#define _ANALYSIS_PARTICLERADIUSDISTRIBUTION_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {
  namespace analysis {
    using namespace iterator;
    /** Class to compute the particle radius distribution. */
    class ParticleRadiusDistribution : public AnalysisBaseTemplate< real > {
    public:
      static void registerPython();

      ParticleRadiusDistribution(shared_ptr< System > system) : AnalysisBaseTemplate< real >(system) {}
      virtual ~ParticleRadiusDistribution() {}

      real computeRaw() {
      int myN, systemN;
      real radsumlocal = 0.0;
      real radsum = 0.0;
      System& system = getSystemRef();
        
      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        radsumlocal += cit->radius() * cit->radius();
      }
          
      myN = system.storage->getNRealParticles();
      
      mpi::all_reduce(*getSystem()->comm, radsumlocal, radsum, std::plus<real>());
      mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());
        
      return 4*radsum/systemN;
      }


      python::list compute() {
        python::list ret;
        real res = computeRaw();
        ret.append(res);
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        real res;
        res = nMeasurements>0 ? newAverage : 0;
        ret.append(res);
        res = nMeasurements>0 ? newVariance : 0;
        ret.append(sqrt(res/(nMeasurements-1)));
        return ret;
      }

      void resetAverage() {
        newAverage   = 0;
        lastAverage  = 0;
        newVariance  = 0;
        lastVariance = 0;
      }

      void updateAverage(real res) {
    	if (nMeasurements > 0) {
    	  if (nMeasurements == 1) {
              newAverage     = res;
              lastAverage    = newAverage;
          } else {
              newAverage   = lastAverage  + (res - lastAverage) / nMeasurements;
              newVariance  = lastVariance + (res - lastAverage) * (res - newAverage);
              lastAverage  = newAverage;
              lastVariance = newVariance;
          }
    	}
        return;
      }
    };
  }
}

#endif
