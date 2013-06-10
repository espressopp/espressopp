// ESPP_CLASS
#ifndef _ANALYSIS_TEMPERATURE_HPP
#define _ANALYSIS_TEMPERATURE_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {
  namespace analysis {
    using namespace iterator;
    /** Class to compute the temperature. */
    class Temperature : public AnalysisBaseTemplate< real > {
    public:
      static void registerPython();

      Temperature(shared_ptr< System > system) : AnalysisBaseTemplate< real >(system) {}
      virtual ~Temperature() {}

      real computeRaw() {
          int myN, systemN;
          real sumT = 0.0;
          real v2sum = 0.0;

          System& system = getSystemRef();

          CellList realCells = system.storage->getRealCells();

          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            Real3D vel = cit->velocity();
            v2sum += cit->mass() * (vel * vel);
          }

          myN = system.storage->getNRealParticles();

          mpi::all_reduce(*getSystem()->comm, v2sum, sumT, std::plus<real>());
          mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());

          return sumT / (3.0 * systemN);
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
