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
#include "GeneralizedLangevinThermostat.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/InterpolationLinear.hpp"
#include "interaction/InterpolationAkima.hpp"
#include "interaction/InterpolationCubic.hpp"


namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;


    GeneralizedLangevinThermostat::GeneralizedLangevinThermostat(shared_ptr<System> system)
    :Extension(system) {

      type = Extension::Thermostat;

      //LOG4ESPP_INFO(theLogger, "GeneralizedLangevinThermostat constructed");

    }

    GeneralizedLangevinThermostat::~GeneralizedLangevinThermostat() {
        disconnect();
    }


    void GeneralizedLangevinThermostat::disconnect() {

        _integrate.disconnect();
        _friction.disconnect();

    }

    void GeneralizedLangevinThermostat::connect() {

        _integrate = integrator->aftIntP.connect(
                boost::bind(&GeneralizedLangevinThermostat::integrate, this));

        _friction = integrator->aftCalcF.connect(
                boost::bind(&GeneralizedLangevinThermostat::friction, this));
    }

    void GeneralizedLangevinThermostat::addCoeffs(int itype, const char* _filename, int type) {
        boost::mpi::communicator world;
        filename = _filename;
        Table table;

        if (itype == 1) { // create a new InterpolationLinear
            table = make_shared <interaction::InterpolationLinear> ();
            table->read(world, _filename);
        }

        else if (itype == 2) { // create a new InterpolationAkima
            table = make_shared <interaction::InterpolationAkima> ();
            table->read(world, _filename);
        }

        else if (itype == 3) { // create a new InterpolationCubic
            table = make_shared <interaction::InterpolationCubic> ();
            table->read(world, _filename);
        }

        coeffs.insert(std::make_pair(type,table));
    }
    
    void GeneralizedLangevinThermostat::integrate()
    {
          //LOG4ESPP_DEBUG(theLogger, "integrate the extended variables");

          System& system = getSystemRef();
          real timestep = integrator->getTimeStep();
          real weight, c, t;
          weight = 0.0;
          c = 0.0;
          t = 0.0;
                  
          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {    
              
              Table table = coeffs.find(cit->getType())->second;
              if (table) { 
                  Particle &vp = *cit;                      
                  weight = vp.lambda();
                  if (weight == 0.0){
                        t = table->getEnergy(weight);
                        vp.extVar() = vp.extVar()*(1.0-timestep/t); 
                  } else if(weight == 1.0){
                        t = table->getEnergy(weight);
                        vp.extVar() = vp.extVar()*(1.0-timestep/t); 
                  } else //if(weight < 1.0 && weight > 0.0)
                  {
                        c = table->getForce(weight);
                        t = table->getEnergy(weight);
                        vp.extVar() = vp.extVar() * (1.0-timestep/t) - (timestep/t) * c *  vp.lambdaDeriv() * vp.lambdaDeriv() * vp.velocity()[0]; 
                  } /*else {
                        std::cout << "ERROR: In Generalized Langevin Thermostat Extension weight is out of [0.0, 1.0] interval: " << weight <<  std::endl;
                        exit(1);
                  }*/                                     
                  
                  //c = table->getForce(weight) * vp.lambdaDeriv * vp.lambdaDeriv;
                  //t = table->getEnergy(weight);
                  //vp.extVar() = vp.extVar()*(1.0-timestep/t) - (timestep/t)*c*vp.velocity()[0]; 
              }
              else{
                  std::cout << "ERROR: Using Generalized Langevin Thermostat Extension without providing table." << std::endl;
                  exit(1);             
              }
              
          }
     }
    
    void GeneralizedLangevinThermostat::friction()
    {
          //LOG4ESPP_DEBUG(theLogger, "apply Generalized Langevin Friction");

          System& system = getSystemRef();

          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
          FixedTupleListAdress::iterator it2;
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {                                      

              Particle &vp = *cit;                      

              it2 = fixedtupleList->find(&vp);
              if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
                    std::vector<Particle*> atList;
                    atList = it2->second;
                    for (std::vector<Particle*>::iterator it3 = atList.begin();
                                         it3 != atList.end(); ++it3) {
                        Particle &at = **it3;
                        at.force()[0] += vp.extVar() * at.mass() / vp.mass();            // GOES BACK IN !!!
                    }
                    //vp.drift() += vp.lambdaDeriv() * fforce;
              }
              else{   // If not, use CG particle itself for calculation.
                         std::cout << "Particle " << vp.id() << " not found in tuples! (Error triggered in GLE extension)" << std::endl;
                         exit(1);
                         return;
              }

          }

    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void GeneralizedLangevinThermostat::registerPython() {


      using namespace espressopp::python;


      class_<GeneralizedLangevinThermostat, shared_ptr<GeneralizedLangevinThermostat>, bases<Extension> >
        ("integrator_GeneralizedLangevinThermostat", init<shared_ptr<System> >())
        .add_property("filename", &GeneralizedLangevinThermostat::getFilename)
        .def("connect", &GeneralizedLangevinThermostat::connect)
        .def("disconnect", &GeneralizedLangevinThermostat::disconnect)
        .def("addCoeffs",  &GeneralizedLangevinThermostat::addCoeffs )// pyAddForce)
        ;

    }

  }
}
