/*
  Copyright (C) 2012,2013,2017,2018
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
#include "mpi.hpp"
#include <mpi4py/mpi4py.h>

namespace espressopp {

  using namespace iterator;
  using std::string;

  namespace integrator {

    ExtPlumed::ExtPlumed(shared_ptr<System> _system, python::object _pyobj, string _plumedfile, string _plumedlog, real _dt):
      Extension(_system),
      plumedfile(_plumedfile),
      plumedlog(_plumedlog),
      dt(_dt),
      step(0),
      nreal(0),
      gatindex(NULL),
      masses(NULL),
      f(NULL),
      pos(NULL),
      charges(NULL)
    {
      System& system = getSystemRef();
      p=new PLMD::Plumed;

      PyObject * pyobj = _pyobj.ptr();
      PyMPICommObject* pyMPIComm = (PyMPICommObject*) pyobj;
      MPI_Comm * comm_p = &pyMPIComm->ob_mpi;
      p->cmd("setMPIComm", comm_p);
      p->cmd("setPlumedDat",plumedfile.c_str());
      p->cmd("setLogFile",plumedlog.c_str());
      p->cmd("setMDEngine","ESPResSo++");
      longint nReal = system.storage->getNRealParticles();
      boost::mpi::all_reduce(*system.comm, nReal, natoms, std::plus<longint>());
      p->cmd("setTimestep",&dt);
      p->cmd("setNatoms",&natoms);
    }

    ExtPlumed::~ExtPlumed() {
      delete [] f;
      delete [] pos;
      delete [] charges;
      delete [] gatindex;
      delete [] masses;
      delete p;
    }

    real ExtPlumed::getBias() {
      return bias;
    }

    void ExtPlumed::setNaturalUnits() {
      p->cmd("setNaturalUnits");
    }

    void ExtPlumed::setTimeUnit(real _factor) {
      p->cmd("setMDTimeUnits",&_factor);
    }

    void ExtPlumed::setLengthUnit(real _factor) {
      p->cmd("setMDLengthUnits",&_factor);
    }

    void ExtPlumed::setEnergyUnit(real _factor) {
      p->cmd("setMDEnergyUnits", &_factor);
    }

    void ExtPlumed::setKbT(real _kbt) {
      p->cmd("setKbT", &_kbt);
    }

    void ExtPlumed::setRealPrecision(int _size) {
      p->cmd("setRealPrecision", &_size);
    }

    void ExtPlumed::setMDChargeUnits(real _factor) {
      p->cmd("setMDChargeUnits", &_factor);
    }

    void ExtPlumed::setMDMassUnits(real _factor) {
      p->cmd("setMDMassUnits", &_factor);
    }

    void ExtPlumed::setRestart(int _res) {
      p->cmd("setRestart", &_res);
    }

    void ExtPlumed::readInputLine(string _input) {
      p->cmd("readInputLine", _input.c_str());
    }

    void ExtPlumed::Init() {
      p->cmd("init");
    }

    void ExtPlumed::disconnect() {
      _runInit.disconnect();
      _aftCalcF.disconnect();
      delete p;
    }

    void ExtPlumed::connect() {
      _runInit = integrator->runInit.connect( boost::bind(&ExtPlumed::setStep, this));
      _aftCalcF = integrator->aftCalcF.connect( boost::bind(&ExtPlumed::updateForces, this));
    }

    void ExtPlumed::setStep() {
      step = integrator->getStep();
    }

    void ExtPlumed::updateForces() {
      int update_gatindex=0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // Try to find out if the domain decomposition has been updated
      if(nreal!=system.storage->getNRealParticles()) {
        if(charges) delete [] charges;
        if(masses) delete [] masses;
        if(gatindex) delete [] gatindex;
        if(pos) delete [] pos;
        if(f) delete [] f;
        nreal = system.storage->getNRealParticles();
        gatindex = new int [nreal];
        masses = new real [nreal];
        charges = new real [nreal];
        pos = new real[nreal*3];
        f = new real[nreal*3];
        update_gatindex=1;
      } else {
        int i = 0;
        for(CellListIterator cit(realCells); !cit.isDone() && i < nreal; ++cit, ++i) {
          if(gatindex[i]!=static_cast<int>(cit->id()-1)) {
            update_gatindex=1;
            break;
          }
        }
      }

      boost::mpi::all_reduce(*system.comm, update_gatindex, std::plus<int>());
      if(update_gatindex) {
        int i=0;
        for(CellListIterator cit(realCells); !cit.isDone() && i < nreal; ++cit, ++i) {
          gatindex[i]=static_cast<int>(cit->id())-1;
          masses[i]=cit->mass();
          charges[i]=cit->q();
          pos[i*3] = cit->position()[0];
          pos[i*3+1] = cit->position()[1];
          pos[i*3+2] = cit->position()[2];
          f[i*3] = cit->force()[0];
          f[i*3+1] = cit->force()[1];
          f[i*3+2] = cit->force()[2];
        }
      } else {
        int i = 0;
        for(CellListIterator cit(realCells); !cit.isDone() && i < nreal; ++cit, ++i) {
          masses[i]=cit->mass();
          charges[i]=cit->q();
          pos[i*3] = cit->position()[0];
          pos[i*3+1] = cit->position()[1];
          pos[i*3+2] = cit->position()[2];
          f[i*3] = cit->force()[0];
          f[i*3+1] = cit->force()[1];
          f[i*3+2] = cit->force()[2];
        }
      }

      p->cmd("setStep",&step);
      p->cmd("setAtomsNlocal",&nreal);
      p->cmd("setAtomsGatindex",&gatindex[0]);
      real box[3][3];

      for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j]=0.0;
      Real3D L = system.bc->getBoxL();
      box[0][0]=L[0];
      box[1][1]=L[1];
      box[2][2]=L[2];

      real virial[3][3];
      for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) virial[i][j]=0.0;

      p->cmd("setPositions", pos);
      p->cmd("setForces", f);
      p->cmd("setBox",&box[0][0]);
      p->cmd("setMasses",masses);
      p->cmd("setCharges",charges);
      p->cmd("setVirial", &virial[0][0]);
      p->cmd("getBias",&bias);

      real pot_energy = 0.;
      const interaction::InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j =0; j < srIL.size(); ++j) {
        pot_energy += srIL[j]->computeEnergy();
      }

      pot_energy /= system.comm->size(); // PLUMED defines PE this way.
      p->cmd("setEnergy", &pot_energy);
      p->cmd("calc");

      // Modifying forces on real particles on each processor.
      int k = 0;
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit, ++k) {
        Real3D F = Real3D(f[k*3], f[k*3+1], f[k*3+2]);
        cit->setF(F);
      }
      step += 1;
    }

    /****************************************************
     ** REGISTRATION WITH PYTHON
     ****************************************************/
    void ExtPlumed::registerPython() {

      using namespace espressopp::python;

      class_<ExtPlumed, shared_ptr<ExtPlumed>, bases<Extension> >

        ("integrator_ExtPlumed", init< shared_ptr< System >, python::object, string, string, real >())
        .def("getBias", &ExtPlumed::getBias)
        .def("setNaturalUnits", &ExtPlumed::setNaturalUnits)
        .def("setTimeUnit", &ExtPlumed::setTimeUnit)
        .def("setEnergyUnit", &ExtPlumed::setEnergyUnit)
        .def("setLengthUnit", &ExtPlumed::setLengthUnit)
        .def("setKbT", &ExtPlumed::setKbT)
        .def("setRealPrecision", &ExtPlumed::setRealPrecision)
        .def("setMDChargeUnit", &ExtPlumed::setMDChargeUnits)
        .def("setMDMassUnit", &ExtPlumed::setMDMassUnits)
        .def("setRestart", &ExtPlumed::setRestart)
        .def("readInputLine", &ExtPlumed::readInputLine)
        .def("Init", &ExtPlumed::Init)
        .def("connect", &ExtPlumed::connect)
        .def("disconnect", &ExtPlumed::disconnect)
        ;
    }
  }
}
