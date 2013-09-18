#include "python.hpp"
#include <cmath>
#include "XPressure.hpp"
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

    // It will calculate the molecular local virial pressure in 'n' layers along X axis, works only for AdResS simulations. 
    python::list XPressure::computeArray(int n) const {
      System& system = getSystemRef();
      //const bc::BC& bc = *system.bc;  // boundary conditions
      system.storage->decompose();

      Real3D Li = system.bc->getBoxL();
      // determine the local volume size
      //real A = Li[1] * Li[2];      
      // n * lX is always Lx
      real lX = Li[0] / (double)n;
      real ThreeV = Li[1] * Li[2] * lX * 3.0;
      int bin = 0;
      real vvlocal[n];
      for(int i=0; i<n; i++) vvlocal[i] = 0.0;
      
      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          
        Particle &vp = *cit;
        FixedTupleListAdress::iterator it2;
        it2 = fixedtupleList->find(&vp);

        if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
              std::vector<Particle*> atList;
              atList = it2->second;

              for (std::vector<Particle*>::iterator it3 = atList.begin();
                                   it3 != atList.end(); ++it3) {
                  
                      Particle &at = **it3;

                      Real3D pos = at.position();
                      real mass = at.mass();
                      Real3D& vel = at.velocity();
                      real vvt = mass * vel[0] * vel[0];

                      if (pos[0] > Li[0])
                      {
                          real pos_wrap = pos[0] - Li[0];
                          bin = floor (pos_wrap / lX);    
                      }         
                      else if (pos[0] < 0.0)
                      {
                          real pos_wrap = pos[0] + Li[0];
                          bin = floor (pos_wrap / lX);    
                      }
                      else
                      {
                          bin = floor (pos[0] / lX);          
                      }
                      
                      vvlocal[bin] += vvt;
                      
                      /*real xminBC = pos[0] - 0.5*lX;
                      real xmaxBC = pos[0] + 0.5*lX;

                      // boundary condition for orthorhombic box
                      bool boundary = false;
                      if(xminBC<0.0){
                        xminBC += Li[0];
                        boundary = true;
                      }
                      else if(xmaxBC>=Li[0]){
                        xmaxBC -= Li[0];
                        boundary = true;
                      }

                      int minpos = (int)( xminBC/lX );
                      int maxpos = (int)( xmaxBC/lX );

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
                      }*/
                      
              }  
        }

        else{   // If not, use CG particle itself for calculation.
              Real3D pos = cit->position();
              real mass = cit->mass();
              Real3D& vel = cit->velocity();
              real vvt = mass * vel[0] * vel[0];
              
              if (pos[0] > Li[0])
              {
                  real pos_wrap = pos[0] - Li[0];
                  bin = floor (pos_wrap / lX);    
              }         
              else if (pos[0] < 0.0)
              {
                  real pos_wrap = pos[0] + Li[0];
                  bin = floor (pos_wrap / lX);    
              }
              else
              {
                  bin = floor (pos[0] / lX);          
              }

              vvlocal[bin] += vvt;
              
              /*real xminBC = pos[0] - 0.5*lX;
              real xmaxBC = pos[0] + 0.5*lX;

              // boundary condition for orthorhombic box
              bool boundary = false;
              if(xminBC<0.0){
                xminBC += Li[0];
                boundary = true;
              }
              else if(xmaxBC>=Li[0]){
                xmaxBC -= Li[0];
                boundary = true;
              }

              int minpos = (int)( xminBC/lX );
              int maxpos = (int)( xmaxBC/lX );

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
              }*/
        }
             
        /*Real3D pos = cit->position();
        real mass = cit->mass();
        Real3D& vel = cit->velocity();
        real vvt = mass * vel[0] * vel[0];
        
        real xminBC = pos[0] - 0.5*lX;
        real xmaxBC = pos[0] + 0.5*lX;
        
        // boundary condition for orthorhombic box
        bool boundary = false;
        if(xminBC<0.0){
          xminBC += Li[0];
          boundary = true;
        }
        else if(xmaxBC>=Li[0]){
          xmaxBC -= Li[0];
          boundary = true;
        }
        
        int minpos = (int)( xminBC/lX );
        int maxpos = (int)( xmaxBC/lX );
        
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
        }*/

      }
      
      real vv[n];
      for ( int i = 0; i < n; ++i)
      {
        vv[i] = 0.0;         
        mpi::all_reduce(*mpiWorld, vvlocal[i], vv[i], std::plus<real>());
      }
      
      // compute the short-range nonbonded contribution
      std::vector<real> w(n);
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialX(w, n);
      }
      
      //python::tuple pijz;
      real wfinal[n];
      python::list XpressureResult;
      for(int i=0; i<n;i++){
        wfinal[i] = w[i];
        vv[i] = vv[i]/ThreeV;       
        XpressureResult.append(vv[i] + wfinal[i]);
      }

      return XpressureResult;
    }
  

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real XPressure::compute() const {
      return -1.0;
    }
    

    using namespace boost::python;

    void XPressure::registerPython() {
      using namespace espresso::python;
      class_<XPressure, bases< Observable > >
        ("analysis_XPressure", init< shared_ptr< System > >())
        .def("compute", &XPressure::computeArray)
      ;
    }
  }
}