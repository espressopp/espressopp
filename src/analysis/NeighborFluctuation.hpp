// ESPP_CLASS
#ifndef _ANALYSIS_NEIGHBORFLUCTUATION_HPP
#define _ANALYSIS_NEIGHBORFLUCTUATION_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"
#include "Real3D.hpp"

#include "storage/DomainDecomposition.hpp"
#include "storage/CellGrid.hpp"

namespace espresso {
  namespace analysis {
    /** Class to get the number of particles in the system. */
    class NeighborFluctuation : public Observable {
    public:
      NeighborFluctuation(shared_ptr< System > system, real _radius) : Observable(system), radius(_radius){
        esutil::Error err(system->comm);
        
        Real3D Li = system->bc->getBoxL();
        
        Int3D cellG = system->storage->getInt3DCellGrid();
        
        real minL = std::min(Li[0]/(real)cellG[0], std::min(Li[1]/(real)cellG[1],Li[2]/(real)cellG[2]));
        if(radius > minL){
          std::stringstream msg;
          msg<<"Error. Radius for checking near neighbors should be smaller then "
                  "minimal cell size.\n";
          msg<<"Otherwise it might count the same particle twice or miss some particles. "
                  "radius="<<radius<<" minCellSize="<<minL;
          err.setException( msg.str() );
        }
      }
      virtual ~NeighborFluctuation() {}
      virtual real compute() const;

      static void registerPython();

    protected:
      real radius;
    };
  }
}

#endif
