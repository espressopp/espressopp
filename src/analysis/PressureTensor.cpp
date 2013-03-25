#include "python.hpp"
#include <cmath>
#include "PressureTensor.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "VerletList.hpp"
#include "storage/NodeGrid.hpp"

#include "Tensor.hpp"

using namespace espresso;
using namespace iterator;
using namespace interaction;

namespace espresso {
  namespace analysis {
    Tensor PressureTensor::computeTensor() const {

      System& system = getSystemRef();

      // determine volume of the box
      Real3D Li = system.bc->getBoxL();
      real V = Li[0] * Li[1] * Li[2];

      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      Tensor vvLocal(0.0);
      Tensor vv;

      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real mass = cit->mass();
        Real3D& vel = cit->velocity();
        vvLocal += mass * Tensor(vel, vel);
      }

      boost::mpi::all_reduce(*mpiWorld, vvLocal, vv, std::plus<Tensor>());

      // compute the short-range nonbonded contribution
      Tensor wij(0.0);
      // virial is already reduced
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialTensor(wij);
      }

      return (vv + wij) / V;
    }

    // it will calculate the pressure in 'n' layers along Z axis
    // the first layer has coordinate Lz/n the last - Lz.
    python::list PressureTensor::computeTensorIKz1(int n, real dz) const {

      System& system = getSystemRef();
      //const bc::BC& bc = *system.bc;  // boundary conditions
      system.storage->decompose();

      Real3D Li = system.bc->getBoxL();
      // determine the local volume size
      real A = Li[0] * Li[1];
      real V = A * (2 * dz);
      
      // n * lZ is always Lz
      real lZ = Li[2] / (double)n;

      Tensor vvlocal[n];
      for(int i=0; i<n; i++) vvlocal[i] = Tensor(0.0);
      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D pos = cit->position();
        real mass = cit->mass();
        Real3D& vel = cit->velocity();
        Tensor vvt = mass * Tensor(vel, vel);
        
        real zminBC = pos[2] - dz;
        real zmaxBC = pos[2] + dz;
        
        // boundary condition for orthorhombic box
        bool boundary = false;
        if(zminBC<0.0){
          zminBC += Li[2];
          boundary = true;
        }
        else if(zmaxBC>=Li[2]){
          zmaxBC -= Li[2];
          boundary = true;
        }
        
        int minpos = (int)( zminBC/lZ );
        int maxpos = (int)( zmaxBC/lZ );
        
        if(boundary){
          
          for(int i = 0; i<=maxpos; i++){
            vvlocal[i] += vvt;
          }
          for(int i = minpos+1; i<n; i++){
            vvlocal[i] += vvt;
          }
        }
        else{
          for(int i = minpos+1; i<=maxpos; i++){
            vvlocal[i] += vvt;
          }
        }

      }
      
      Tensor vv[n];
      mpi::all_reduce(*mpiWorld, vvlocal, n, vv, std::plus<Tensor>());
      
      // compute the short-range nonbonded contribution
      Tensor w[n];
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialTensor(w, n);
      }
      
      //python::tuple pijz;
      python::list pijarr;
      for(int i=0; i<n;i++){
        vv[i] = vv[i]/V;
        w[i] = w[i] / A;
        
        pijarr.append(vv[i] + w[i]);
      }

      return pijarr;
    }

    // calculates the pressure for layer with z coordinate
    // kinetic part: particles in [z-dz, z+dz] slab are taken into account
    Tensor PressureTensor::computeTensorIKz2(real z, real dz) const {
      System& system = getSystemRef();
      Real3D Li = system.bc->getBoxL();

      // determine the local volume size for kinetic part
      real A = Li[0] * Li[1];
      real V = A * (2 * dz);

      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      Tensor vvlocal(0.0);
      CellList realCells = system.storage->getRealCells();
      
      real zmin = z - dz;
      real zmax = z + dz;
      
      real zminBC = zmin;
      real zmaxBC = zmax;
      
      // boundary condition for orthorhombic box
      bool condition2 = false;
      if(zminBC<0.0){
        zminBC += Li[2];
        condition2 = true;
      }
      else if(zmaxBC>=Li[2]){
        zmaxBC -= Li[2];
        condition2 = true;
      }
      
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D pos = cit->position();
        
        bool condition;
        if(condition2){
          condition = (pos[2]>=zminBC || pos[2]<zmaxBC);
        }
        else{
          condition = (pos[2]>=zminBC && pos[2]<zmaxBC);
        }
        
        if( condition ){
          real mass = cit->mass();
          Real3D& vel = cit->velocity();
          vvlocal += mass * Tensor(vel, vel);
        }
      }
      Tensor vv(0.0);
      boost::mpi::all_reduce(*mpiWorld, vvlocal, vv, std::plus<Tensor>());
      
      vv = vv / V;

      // compute the short-range nonbonded contribution
      Tensor w(0.0);
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialTensor(w, z);
      }
      
      w = w / A;
      
      return ( vv + w );
    }

    real PressureTensor::compute() const {
      return -1.0;
    }

    // TODO it is fast solution. one should think about the overloading
    
    using namespace boost::python;

    void PressureTensor::registerPython() {
      using namespace espresso::python;
      class_<PressureTensor, bases< Observable > >
        ("analysis_PressureTensor", init< shared_ptr< System > >())
        .def("compute1", &PressureTensor::computeTensor)
        .def("compute2", &PressureTensor::computeTensorIKz1)
        .def("compute3", &PressureTensor::computeTensorIKz2)
      ;
    }
  }
}
