/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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
#include "FreeEnergyCompensation.hpp"
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

    LOG4ESPP_LOGGER(FreeEnergyCompensation::theLogger, "FreeEnergyCompensation");

    FreeEnergyCompensation::FreeEnergyCompensation(shared_ptr<System> system, bool _sphereAdr, int _ntrotter, bool _slow)
    :Extension(system), sphereAdr(_sphereAdr), ntrotter(_ntrotter), slow(_slow) {

        type = Extension::FreeEnergyCompensation;

        center = (0.0,0.0,0.0);

        LOG4ESPP_INFO(theLogger, "FreeEnergyCompensation constructed");
    }

    FreeEnergyCompensation::~FreeEnergyCompensation() {}



    void FreeEnergyCompensation::connect(){
      if(slow){
        _applyForce = integrator->aftCalcSlow.connect(
            boost::bind(&FreeEnergyCompensation::applyForce, this));
      }
      else{
        _applyForce = integrator->aftCalcF.connect(
            boost::bind(&FreeEnergyCompensation::applyForce, this));
      }
    }

    void FreeEnergyCompensation::disconnect(){
        _applyForce.disconnect();
    }


    void FreeEnergyCompensation::addForce(int itype, const char* _filename, int type) {
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

        forces.insert(std::make_pair(type,table));
    }


    void FreeEnergyCompensation::applyForce() {
          LOG4ESPP_DEBUG(theLogger, "apply Free Energy Compensation force");

          System& system = getSystemRef();

          std::unordered_map<int, Table>::iterator tableIt;
          Table table;
          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
          FixedTupleListAdress::iterator it2;
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
              tableIt = forces.find(cit->getType());
              if (tableIt != forces.end()) { //because there may be CG particles to which Free Energy Compensation is not applied
                table = tableIt->second;
              }

              if (table) {

                  Particle &vp = *cit;
                  real weight = vp.lambda();

                  if (weight != 1.0 && weight != 0.0){
                          real fforce = table->getForce(weight);

                          Real3D dist3D(0.0, 0.0, 0.0);
                          if(sphereAdr){
                            dist3D = cit->getPos() - center;
                            real absdist3D = sqrt(dist3D.sqr());
                            dist3D *= fforce;
                            dist3D /= absdist3D;
                          }
                          else{
                            real dist = vp.position()[0]-center[0];
                            if ( dist >= 0.0 ) {dist = 1.0;}
                            else {dist = -1.0;}
                            fforce *=dist;
                          }

                          it2 = fixedtupleList->find(&vp);
                          if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
                                std::vector<Particle*> atList;
                                atList = it2->second;
                                for (std::vector<Particle*>::iterator it3 = atList.begin();
                                                     it3 != atList.end(); ++it3) {
                                    Particle &at = **it3;

                                    if(sphereAdr){
                                      at.force() += vp.lambdaDeriv() * at.mass() * dist3D / (vp.mass() * ntrotter);
                                    }
                                    else{
                                      at.force()[0] += vp.lambdaDeriv() * at.mass() * fforce / (vp.mass() * ntrotter);
                                    }

                                }
                          }
                          else{   // If not, use CG particle itself for calculation.
                                     std::cout << "Particle " << vp.id() << " not found in tuples!" << std::endl << "It's unclear how FEC work when combining particles, which do change resolution with particles that don't." << std::endl;
                                     exit(1);
                                     return;
                          }
                  }
              }
          }
    }

    real FreeEnergyCompensation::computeCompEnergy() {
          LOG4ESPP_DEBUG(theLogger, "compute Free Energy Compensation Energies");

          real CompEnergy = 0.0;
          real CompEnergySum = 0.0;
          System& system = getSystemRef();

          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

              Table table = forces.find(cit->getType())->second;
              if (table) {
                  Particle &vp = *cit;
                  real weight = vp.lambda();
                  CompEnergy += table->getEnergy(weight);
              }
              else{
                  std::cout << "ERROR: Using FEC Extension without providing table." << std::endl;
                  exit(1);
                  return 0.0;
              }

          }
          mpi::all_reduce(*getSystem()->comm, CompEnergy, CompEnergySum, std::plus<real>());
          return CompEnergySum;
    }

    void FreeEnergyCompensation::setCenter(real x, real y, real z){
            center = Real3D(x, y, z);
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FreeEnergyCompensation::registerPython() {

      using namespace espressopp::python;

      class_<FreeEnergyCompensation, shared_ptr<FreeEnergyCompensation>, bases<Extension> >
        ("integrator_FreeEnergyCompensation", init< shared_ptr<System>, bool, int, bool >())
        .add_property("filename", &FreeEnergyCompensation::getFilename)
        .def("connect", &FreeEnergyCompensation::connect)
        .def("disconnect", &FreeEnergyCompensation::disconnect)
        .def("setCenter", &FreeEnergyCompensation::setCenter) // pySetCenter)
        .def("addForce",  &FreeEnergyCompensation::addForce )// pyAddForce)
        .def("computeCompEnergy", &FreeEnergyCompensation::computeCompEnergy)
        ;
    }

  }
}
