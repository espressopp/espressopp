/*
 Copyright (C) 2012-2016
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

#include "LBOutputScreen.hpp"
#include "storage/Storage.hpp"

namespace espressopp {
   namespace analysis {
      //  LOG4ESPP_LOGGER(LBOutputScreen::theLogger, "LBOutputScreen");
      LBOutputScreen::LBOutputScreen(shared_ptr<System> system,
                                     shared_ptr< integrator::LatticeBoltzmann >latticeboltzmann)
      : LBOutput(system, latticeboltzmann) {}

      void LBOutputScreen::writeOutput() {
         /* It is a not straightforward way. The idea is that the Output onto
          the Screen takes place when the step is already finished, i.e. at
          (t+ dt). However, the LB couples to the integrator signal befIntV,
          i.e. at (t+1/2dt). Therefore, to check momentum conservation for LB+MD,
          it is necessary to go forward in time by 1/2dt in MD-part. This is
          done by employing the same forces as for integrate1() but withont
          actually altering velocities of MD-particles. */

         bool _restart = latticeboltzmann->doRestart();
         real _timestep = latticeboltzmann->getCopyTimestep();

         Real3D _myCMVel = Real3D(0.);               // cm velocity (current CPU)
         Real3D _totCMVel = Real3D(0.);            // cm velocity (all CPUs)

         System& system = getSystemRef();
         CellList realCells = system.storage->getRealCells();

         // loop over all particles in the curr CPU //
         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            Real3D& vel = cit->velocity();
            Real3D& force = cit->force();
            _myCMVel += vel + 0.5 * _timestep * force;
         }

         if (_restart) {
            // collect cm vels and share the result between CPUs
            boost::mpi::all_reduce(*system.comm,
                                   _myCMVel, _totCMVel, std::plus<Real3D>());
            findLBMom(1);

            // correct for difference between LB and MD momenta
            Real3D corrForCMVel = getLBMom() + _totCMVel;
            Real3D corrPerPart = corrForCMVel / latticeboltzmann->getTotNPart();
            latticeboltzmann->galileanTransf(corrPerPart);

            // set flag after restart
            latticeboltzmann->setDoRestart(false);
         } else {
            // collect cm vels from CPUs
            boost::mpi::reduce(*system.comm,
                               _myCMVel, _totCMVel, std::plus<Real3D>(), 0);
            findLBMom(0);
         }

         // output MD momentum and statistics of the run //
         if (system.comm->rank() == 0) {
            Int3D _Ni = latticeboltzmann->getNi();
            long int latSize = _Ni[0] * _Ni[1] * _Ni[2];
            long int _step = latticeboltzmann->getStepNum();
            real _invNSteps = 1. / latticeboltzmann->getNSteps();
            long int madeLBSteps = (long int)((_step-getOldStepNum()) * _invNSteps);

            // output MD momentum //
            setMDMom(_totCMVel); // need this for testsuite
            printf("MD mom at (+1/2dt): %17.12f %17.12f %17.12f \n",
                   _totCMVel[0], _totCMVel[1], _totCMVel[2]);

            if (_step != 0) {
               // calculate time performance
               setLBTimerNew( timeLBtoMD.getElapsedTime() );
               real timelb = getLBTimerNew() - getLBTimerOld();
               setLBTimerOld( getLBTimerNew() );

               // output LB statistics
               printf ("%ld LB steps on %ld nodes done in %f sec, relative MLUPS: %f \n",
                       madeLBSteps, latSize, timelb,
                       madeLBSteps * latSize * 1e-6 / timelb);
               printf ("Speed in LJ units: %f tau / sec \n",
                       (_step - getOldStepNum()) * _timestep / timelb);

            } else {
               // timers //
               timeLBtoMD.reset();
               setLBTimerOld(0.);
            }

            setOldStepNum(_step);

            printf("------------------------------------");
            printf("-------------------------------------\n");
         }
      }

      void LBOutputScreen::findLBMom (int _mode) {
         Real3D result = Real3D(0.);
         Real3D _myMom = Real3D(0.);
         Int3D _Ni = latticeboltzmann->getMyNi();
         int _offset = latticeboltzmann->getHaloSkin();

         for (int i = _offset; i < _Ni[0] - _offset; i++) {
            for (int j = _offset; j < _Ni[1] - _offset; j++) {
               for (int k = _offset; k < _Ni[2] - _offset; k++) {
                  _myMom += Real3D( latticeboltzmann->getLBMom(Int3D(i,j,k), 1),
                                    latticeboltzmann->getLBMom(Int3D(i,j,k), 2),
                                    latticeboltzmann->getLBMom(Int3D(i,j,k), 3) );
               }
            }
         }
         _myMom *= (latticeboltzmann->convTimeMDtoLB() / latticeboltzmann->getA());

         System& system = getSystemRef();

         if (_mode == 0) {
            // collect momenta from CPUs
            boost::mpi::reduce(*system.comm,
                               _myMom, result, std::plus<Real3D>(), 0);
         } else {
            // collect momenta and share the result between CPUs
            boost::mpi::all_reduce(*system.comm,
                                   _myMom, result, std::plus<Real3D>());
            setLBMom(result);
         }

         // output LB momentum //
         if (system.comm->rank() == 0) {
            setLBMom(result); // need this for testsuite
            printf ("after step %d:\n", latticeboltzmann->getStepNum());
            printf ("LB mom in LJ units: %17.12f %17.12f %17.12f \n",
                    result[0], result[1], result[2]);
         }
      }

      void LBOutputScreen::setLBMom (Real3D _lbMom) { lbMom = _lbMom; }
      Real3D LBOutputScreen::getLBMom () {return lbMom;}

      void LBOutputScreen::setMDMom (Real3D _mdMom) { mdMom = _mdMom; }
      Real3D LBOutputScreen::getMDMom () {return mdMom;}

      void LBOutputScreen::setLBTimerOld (real _lbTimerOld) {
         lbTimerOld = _lbTimerOld;}
      real LBOutputScreen::getLBTimerOld() {return lbTimerOld;}

      void LBOutputScreen::setLBTimerNew(real _lbTimerNew) {
         lbTimerNew = _lbTimerNew;}
      real LBOutputScreen::getLBTimerNew() {return lbTimerNew;}

      void LBOutputScreen::setOldStepNum(long int _oldStepNum) {
         oldStepNum = _oldStepNum;}
      long int LBOutputScreen::getOldStepNum() {return oldStepNum;}

      void LBOutputScreen::registerPython() {
         using namespace espressopp::python;

         class_<LBOutputScreen, bases< LBOutput > >
         ("analysis_LBOutput_Screen", init< shared_ptr< System >,
          shared_ptr< integrator::LatticeBoltzmann > >())

         .def("writeOutput", &LBOutputScreen::writeOutput)
         .def("getLBMom", &LBOutputScreen::getLBMom)
         .def("getMDMom", &LBOutputScreen::getMDMom)
         ;
      }
   }
}

