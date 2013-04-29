#include "python.hpp"
#include <cmath>
#include "Temperature.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "FixedTupleList.hpp"

using namespace espresso;
using namespace iterator;

namespace espresso {
  namespace analysis {
    real Temperature::compute() const {

      int myN, systemN;
      real sumT = 0.0;
      real v2sum = 0.0;
      System& system = getSystemRef();
        
      if (system.storage->getFixedTuples()){  // AdResS - hence, need to distinguish between CG and AT particles.     
          int count = 0;
          shared_ptr<FixedTupleList> fixedtupleList=system.storage->getFixedTuples();
          CellList realCells = system.storage->getRealCells();
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {  // Iterate over all (CG) particles.              
            Particle &vp = *cit;
            FixedTupleList::iterator it2;
            it2 = fixedtupleList->find(&vp);
            
            if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
                  std::vector<Particle*> atList;
                  atList = it2->second;
                  for (std::vector<Particle*>::iterator it3 = atList.begin();
                                       it3 != atList.end(); ++it3) {
                      Particle &at = **it3;
                      Real3D vel = at.velocity();
                      v2sum += at.mass() * (vel * vel);
                      count += 1;
                  }  
            }
            
            else{   // If not, use CG particle itself for calculation.
                  Real3D vel = cit->velocity();
                  v2sum += cit->mass() * (vel * vel);
                  count += 1;
            }
            
          }
          
          myN = count;
          //std::cout << "COUNT " << count << std::endl;
      }
      
      else{  // No AdResS - just iterate over all particles          
          CellList realCells = system.storage->getRealCells();
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            Real3D vel = cit->velocity();
            v2sum += cit->mass() * (vel * vel);
          }
          
          myN = system.storage->getNRealParticles();
      }
      
      mpi::all_reduce(*getSystem()->comm, v2sum, sumT, std::plus<real>());
      mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());

      return sumT / (3.0 * systemN); 
      
    }
    
    void Temperature::registerPython() {
      using namespace espresso::python;
      class_<Temperature, bases< Observable > >
        ("analysis_Temperature", init< shared_ptr< System > >())
      ;
    }
  }
}
