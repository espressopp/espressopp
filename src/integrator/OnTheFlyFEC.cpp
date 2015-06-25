/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"
#include "OnTheFlyFEC.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/InterpolationLinear.hpp"
#include "interaction/InterpolationAkima.hpp"
#include "interaction/InterpolationCubic.hpp"


namespace espressopp {
  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(OnTheFlyFEC::theLogger, "OnTheFlyFEC");

    OnTheFlyFEC::OnTheFlyFEC(shared_ptr<System> system)
    :Extension(system) {

        type = Extension::FreeEnergyCompensation;
        
        bins = 0;
        gap = 0;
        steps = 0;
        center = (0.0,0.0,0.0);
        
        counter = 0;
        gapcounter = 0;
        
        EnergyDiff = 0;
        NumbersAtoms = 0;

        LOG4ESPP_INFO(theLogger, "OnTheFlyFEC constructed");
    }

    OnTheFlyFEC::~OnTheFlyFEC() {
            disconnect();
    }

    // Connect & Disconnect
    void OnTheFlyFEC::connect(){
        _gatherStats = integrator->aftIntV.connect(
            boost::bind(&OnTheFlyFEC::gatherStats, this));
    }

    void OnTheFlyFEC::disconnect(){
         delete[] EnergyDiff;
         EnergyDiff = 0;
         delete[] NumbersAtoms;
         NumbersAtoms = 0;
        _gatherStats.disconnect();
    }

    // Setter & Getter
    void OnTheFlyFEC::setBins(int _bins)
    {
      bins = _bins;
    }
    int OnTheFlyFEC::getBins()
    {
      return bins;
    }

    void OnTheFlyFEC::setGap(int _gap)
    {
      gap = _gap;
    }
    int OnTheFlyFEC::getGap()
    {
      return gap;
    }

    void OnTheFlyFEC::setSteps(int _steps)
    {
      steps = _steps;
    }
    int OnTheFlyFEC::getSteps()
    {
      return steps;
    }
    
    // Make the Vectors
    void OnTheFlyFEC::makeArrays()
    {
      NumbersAtoms = new int[bins];
      for(int i=0;i<bins;i++) NumbersAtoms[i]=0;
    
      EnergyDiff = new real[bins];
      for(int i=0;i<bins;i++) EnergyDiff[i]=0;    
        
      /*for(int i = 0; i < bins; i++)
      { 
        NumbersAtoms.push_back(0); 
        EnergyDiff.push_back(0.0);
      }*/ 
    }
    
    // Gather Statistics for "steps" steps    
    void OnTheFlyFEC::gatherStats() {
          LOG4ESPP_DEBUG(theLogger, "gather Stats for OnTheFlyFEC");
          
          System& system = getSystemRef();

          gapcounter += 1;

          // Only go in after gap steps
					// Nikita: just changed the '=' to '==' as I think it was initially meant
          if (gapcounter == gap){
          
              // Only go in if counter below steps
              if (counter < (steps/gap)){
                  
                  // iterate over CG particles
                  CellList cells = system.storage->getRealCells();
                  shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
                  FixedTupleListAdress::iterator it2;
                  for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

                      Particle &vp = *cit;
                      real weight = vp.lambda();  

                      if (weight < 0.9999999 && weight > 0.0000001){
                          //std::cout << "weight: " << weight << "\n";
                          //std::cout << "floor(weight*bins): " << floor(weight*bins) << "\n";

                          real driftterm = vp.drift();

                          real dist = vp.position()[0]-center[0];
                          if ( dist < 0.0 ) {driftterm = (-1.0)*driftterm ;}     

                          NumbersAtoms[int(floor(weight*bins))] += 1;
                          EnergyDiff[int(floor(weight*bins))] += driftterm;                                    
                      }
                  }
                  counter += 1;
                  
              }
              
              gapcounter = 0;
              
          }
          
    }
    
    // Reset the counter
    void OnTheFlyFEC::resetCounter()
    {
        counter = 0;
    }
   
    // Write the current data
    python::list OnTheFlyFEC::writeFEC() {
          LOG4ESPP_DEBUG(theLogger, "write OnTheFlyFEC");
          
      //std::vector<int> NumbersAtomsTotal;
      //std::vector<real> EnergyDiffTotal;

      int * NumbersAtomsTotal = 0;
      NumbersAtomsTotal = new int[bins];
      for(int i=0;i<bins;i++) NumbersAtomsTotal[i]=0;
      
      real * EnergyDiffTotal = 0;
      EnergyDiffTotal = new real[bins];
      for(int i=0;i<bins;i++) EnergyDiffTotal[i]=0.0;      
      
      /*for(int i = 0; i < bins; i++)
      { 
        NumbersAtomsTotal.push_back(0); 
        EnergyDiffTotal.push_back(0.0);
      }*/ 

      /*for(int i=0; i < bins; i++){
        std::cout << "i: " << i << "\n";
        std::cout << "EnergyDiffTotal.at(i): " << &(EnergyDiffTotal.at(i)) << "\n";
        std::cout << "NumbersAtomsTotal.at(i): " << &(NumbersAtomsTotal.at(i)) << "\n";
        std::cout << "EnergyDiffTotal[i]: " << &(EnergyDiffTotal[i]) << "\n";
        std::cout << "NumbersAtomsTotal[i]: " << &(NumbersAtomsTotal[i]) << "\n";
      }*/
      //real * totHistogram = 0;
      //totHistogram = new real[splitN];
      //for(int i=0;i<splitN;i++) totHistogram[i]=0.0;
     
      //boost::mpi::reduce(*mpiWorld, (real*)&EnergyDiff, bins, (real*)&EnergyDiffTotal, std::plus<real>(),0);
      //boost::mpi::reduce(*mpiWorld, (int*)&NumbersAtoms, bins, (int*)&NumbersAtomsTotal, std::plus<int>(),0);
      boost::mpi::all_reduce(*mpiWorld, EnergyDiff, bins, EnergyDiffTotal, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, NumbersAtoms, bins, NumbersAtomsTotal, std::plus<int>());

      // normalizing
      //int nconfigs = 1; //config - 1
      //real rho = (real)total_num_part / splitN;
      
      //for(int i=0; i < splitN; i++){
      //  totHistogram[i] /= rho;
      //}

      /*std::cout << "bins: " << bins << "\n";
      for(int i=0; i < bins; i++){
        std::cout << "i: " << i << "\n";
        std::cout << "EnergyDiff.at(i): " << EnergyDiff.at(i) << "\n";
        std::cout << "NumbersAtoms.at(i): " << NumbersAtoms.at(i) << "\n";
      }
      
      std::cout << "bins: " << bins << "\n";
      for(int i=0; i < bins; i++){
        std::cout << "i: " << i << "\n";
        std::cout << "EnergyDiffTotal.at(i): " << &(EnergyDiffTotal.at(i)) << "\n";
        std::cout << "NumbersAtomsTotal.at(i): " << &(NumbersAtomsTotal.at(i)) << "\n";
        std::cout << "EnergyDiffTotal[i]: " << &(EnergyDiffTotal[i]) << "\n";
        std::cout << "NumbersAtomsTotal[i]: " << &(NumbersAtomsTotal[i]) << "\n";
      }*/
      
      python::list pyli;
      for(int i=0; i < bins; i++){
        //std::cout << "Loop iter: " << i << "\n";
          if (NumbersAtomsTotal[i] > 0){
                real x = EnergyDiffTotal[i]/NumbersAtomsTotal[i]; 
                pyli.append(x);
                //std::cout << "EnergyDiffTotal[" << i << "] = " << EnergyDiffTotal[i] << "\n";
                //std::cout << "NumbersAtomsTotal[" << i << "] = " << NumbersAtomsTotal[i] << "\n";
                //std::cout << "EnergyDiffTotal/NumbersAtomsTotal[" << i << "] = " << x << "\n";
          }
          else{
                pyli.append(0.0);
                std::cout << "Warning: No sampling for OnTheFlyFEC in bin " << i << " done. Appended 0.0\n"; 
          }
      }
      //std::cout << "Last Loop done " << "\n";
      //EnergyDiffTotal.clear();
      //NumbersAtomsTotal.clear();
      //std::cout << "Maps cleared " << "\n";
      //for(int i=0; i < bins; i++){
        //std::cout << "pyli[" << i << "] = " << pyli[i] << "\n";
        //  std::cout << pyli << "\n";
      //}
      delete[] EnergyDiffTotal;
      EnergyDiffTotal = 0;
      delete[] NumbersAtomsTotal;
      NumbersAtomsTotal = 0;

      return pyli;    
    }
    
    void OnTheFlyFEC::setCenter(real x, real y, real z){
          center = Real3D(x, y, z);
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void OnTheFlyFEC::registerPython() {

      using namespace espressopp::python;


      class_<OnTheFlyFEC, shared_ptr<OnTheFlyFEC>, bases<Extension> >
        ("integrator_OnTheFlyFEC", init< shared_ptr<System> >())
        .add_property("bins", &OnTheFlyFEC::getBins, &OnTheFlyFEC::setBins)
        .add_property("steps", &OnTheFlyFEC::getSteps, &OnTheFlyFEC::setSteps)
        .add_property("gap", &OnTheFlyFEC::getGap, &OnTheFlyFEC::setGap)
        .def("connect", &OnTheFlyFEC::connect)
        .def("disconnect", &OnTheFlyFEC::disconnect)
        .def("makeArrays", &OnTheFlyFEC::makeArrays) // pyMakeVecs)
        .def("writeFEC",  &OnTheFlyFEC::writeFEC )// pyWriteFec)
        .def("resetCounter", &OnTheFlyFEC::resetCounter)
        .def("getBins", &OnTheFlyFEC::getBins)
        .def("getSteps", &OnTheFlyFEC::getSteps)
        .def("getGap", &OnTheFlyFEC::getGap)
        .def("setCenter", &OnTheFlyFEC::setCenter) // pySetCenter)
        ;
    }

  }
}
