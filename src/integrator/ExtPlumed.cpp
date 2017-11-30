/*
  Copyright (C) 2012,2013,2017
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

#include <string>
#include "python.hpp"
#include "types.hpp"
#include "Real3D.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "ExtPlumed.hpp"
// #include "mpi.h"
/*
Todo: units. Not just natural units.
Todo: particle groups. ExtPlumed::applyForceToGroup()
*/
namespace espressopp {

  using namespace iterator;
  using namespace std;

  namespace integrator {

    LOG4ESPP_LOGGER(ExtPlumed::theLogger, "ExtPlumed");

    ExtPlumed::ExtPlumed(shared_ptr<System> _system, string filename, string plumedlog, string units)
      : Extension(_system),
	plumedfile(filename),
        logfile(plumedlog),
        units(units),
        // pe(nullptr),
        nlocal(0),
        natoms(0),
        charged(false),
        gatindex(nullptr),
        masses(nullptr),
        forces(nullptr),
        pos(nullptr),
        charges(nullptr)
    {
      System& system = getSystemRef();
      p=new PLMD::Plumed();

      // p->cmd("setMPIComm",*system.comm);

      // Todo: support more units
      if (units == "Natural") p->cmd("setNaturalUnits");

      p->cmd("setPlumedDat",plumedfile.c_str());
      p->cmd("setLogFile",logfile.c_str());
      p->cmd("setMDEngine","ESPRESSO++");

      longint nReal = system.storage->getNRealParticles();
      // boost::mpi::all_reduce(*system.comm, nReal, natoms, plus<longint>());
      natoms = nReal;
      p->cmd("setNatoms",&natoms);  // check type
      dt = integrator->getTimeStep();
      p->cmd("setTimestep",&dt); // check type
      p->cmd("init");

      // auto IL = system.getInteraction(0);
      // pe = new analysis::PotentialEnergy(system, system.getInteraction(0));
      nprocs = system.comm->size();
      // peTotal = pe->compute_real() / nprocs; // pe is defined as per proc in plumed.
      // p->cmd("setEnergy", &peTotal);
      p->cmd("setAtomsNlocal",&nlocal);
    }

        real ExtPlumed::getBias(){
          return bias;
        }

    void ExtPlumed::disconnect(){
      _aftCalcF.disconnect();
      // _aftIntV.disconnect();
    }

    void ExtPlumed::connect(){
      // connection to initialisation
      _aftCalcF  = integrator->aftCalcF.connect( boost::bind(&ExtPlumed::applyForceToAll, this));
      // _aftIntV = integrator->aftIntV.connect( boost::bind(&ExtPlumed::computePe, this));
    }

    void ExtPlumed::applyForceToAll() {
      int update_gatindex=0;
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      // Try to find out if the domain decomposition has been updated:
      if(nlocal!=system.storage->getNLocalParticles()) {
        if(charges) delete [] charges;
        if(masses) delete [] masses;
        if(gatindex) delete [] gatindex;
        if(forces) delete [] forces;
        if(pos) delete [] pos;
        nlocal = system.storage->getNLocalParticles();
        gatindex = new size_t[nlocal];
        masses = new real[nlocal];
        forces = new real[nlocal*3];
        pos = new real[nlocal*3];
        if (charged) charges = new real[nlocal];
        update_gatindex=1;
      } else {
        longint pcount = 0;
        for(CellListIterator cit(localCells); !cit.isDone() && pcount < nlocal; ++cit, ++pcount) {
          if(gatindex[pcount]!=cit->id()-1) {
            update_gatindex=1;
            break;
          }
          // Check if the number of local particles of plumed is equal to the
          // total number of particles stored in all cells.
          assert(pcount==nlocal);
        }
      }

      // boost::mpi::all_reduce(*system.comm, update_gatindex, plus<int>());
      if(update_gatindex) {
        int i=0;
        for(CellListIterator cit(localCells); !cit.isDone() && i < nlocal; ++cit) {
          gatindex[i]=cit->id()-1;
          masses[i]=cit->mass();
          if (charged) charges[i]=cit->q(); // maybe aslo charges.
          pos[3*i] = cit->position()[0];
          pos[3*i+1] = cit->position()[1];
          pos[3*i+2] = cit->position()[2];
          forces[3*i] = cit->force()[0];
          forces[3*i+1] = cit->force()[1];
          forces[3*i+2] = cit->force()[2];
          ++i;
        }
      } else {
        int i = 0;
        for(CellListIterator cit(localCells); !cit.isDone() && i < nlocal; ++cit) {
          pos[3*i] = cit->position()[0];
          pos[3*i+1] = cit->position()[1];
          pos[3*i+2] = cit->position()[2];
          forces[3*i] = cit->force()[0];
          forces[3*i+1] = cit->force()[1];
          forces[3*i+2] = cit->force()[2];
          ++i;
        }
      }
      p->cmd("setAtomsGatindex",&gatindex[0]);

      real box[3][3];
      for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j]=0.0;
      Real3D L = system.bc->getBoxL();
      box[0][0]=L[0];
      box[1][1]=L[1];
      box[2][2]=L[2];

      // local variable with timestep:
      // step: long long
      int step=integrator->getStep();
      p->cmd("setStep",&step);
      p->cmd("setPositions",&pos[0]);
      p->cmd("setBox",&box[0][0]);
      p->cmd("setForces",&forces[0]);
      p->cmd("setMasses",&masses[0]);
      p->cmd("setCharges",&charges[0]);
      p->cmd("getBias",&bias);
      p->cmd("calc");

      int k =0;
      for(CellListIterator cit(localCells); !cit.isDone() && k < nlocal; ++cit) {
        Real3D f = Real3D(forces[3*k], forces[3*k+1], forces[3*k+2]);
        cit->setF(f);
        ++k;
      }
    }

    // void ExtPlumed::computePe() {
    //   peTotal = pe->compute_real() / nprocs;
    //   return;
    // }

    // ExtPlumed::~ExtPlumed() {
    //   if(pe != nullptr) {
    //       delete pe;
    //   }
    // }

    /****************************************************
     ** REGISTRATION WITH PYTHON
     ****************************************************/
    void ExtPlumed::registerPython() {

      using namespace espressopp::python;

      class_<ExtPlumed, shared_ptr<ExtPlumed>, bases<Extension> >

        ("integrator_ExtPlumed", init< shared_ptr< System>, string, string, string >())
        .def("connect", &ExtPlumed::connect)
        .def("disconnect", &ExtPlumed::disconnect)
        ;
    }
  }
}
