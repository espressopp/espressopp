/*
 Copyright (C) 2012-2016 Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
 
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
//    LOG4ESPP_LOGGER(LBOutputScreen::theLogger, "LBOutputScreen");
    LBOutputScreen::LBOutputScreen(shared_ptr<System> system,
                                       shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
    : LBOutput(system, latticeboltzmann) {}

    void LBOutputScreen::writeOutput()
    {
			/* GET VELOCITIES FROM PREVIOUS HALF-TIMESTEP */
			/* be careful, it is a bit "dirty" way. The idea is that the Output onto the Screen
			 takes place when the timestep is already finished, i.e. at (t+dt). However, the LB
			 algorithm connects to the integrator at (t+1/2dt). Therefore, to check momentum
			 conservation for LB+MD, it is necessary to go forward in time by 1/2dt in MD-part.
			 This is done by using the same forces as for integrate1() but withont actually 
			 altering velocities of MD-particles. */
			real _timestep = latticeboltzmann->getCopyTimestep();
			
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
			
			Real3D _velCM = Real3D(0.,0.,0.);
			Real3D _totVelCM = Real3D(0.,0.,0.);
			
			// loop over all particles in the curr CPU and interpolate CM's vel to moment t+dt
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				Real3D& vel = cit->velocity();
				Real3D& force = cit->force();
				_velCM += vel + 0.5 * _timestep * force;
			}
			
			boost::mpi::reduce(*getSystem()->comm, _velCM, _totVelCM, std::plus<Real3D>(), 0);
			
			findLBMom();
	
			if (getSystem()->comm->rank() == 0) {
				Int3D _Ni = latticeboltzmann->getNi();
				long int lbVolume = _Ni[0]*_Ni[1]*_Ni[2];
				
				printf("_velCM : cmV(t+ 1/2dt) of LJ system is     %18.14f %18.14f %18.14f \n",
							 _totVelCM[0], _totVelCM[1], _totVelCM[2]);
				
				long int _step = latticeboltzmann->getStepNum();
            real _invNSteps = 1. / latticeboltzmann->getNSteps();
				
				if (_step != 0) {
					// calculate time performance
					real timelb;
					setLBTimerNew(timeLBtoMD.getElapsedTime());
					timelb = getLBTimerNew() - getLBTimerOld();
					setLBTimerOld(getLBTimerNew());
					printf ("_step is %ld, getOldStepNum() is %ld, lbVolume is %ld\n", _step, getOldStepNum(), lbVolume);
					printf ("time spent on %ld LB steps is %f sec, relative MLUPS: %f \n",
									(long int)((_step-getOldStepNum())*_invNSteps), timelb,
									(_step-getOldStepNum())*_invNSteps*lbVolume*1e-6 / timelb);
					printf ("Speed of simulation (in LJ units): %f tau per second\n",
									(_step-getOldStepNum())*_timestep / timelb);

				} else {
					timeLBtoMD.reset();
					setLBTimerOld(0.);
					printf ("Initialisation of the LB-MD system has finished\n");
				}

				setOldStepNum(_step);

				printf("-------------------------------------------------------------------\n");
			}
    }
		
		void LBOutputScreen::findLBMom () {
			real _fi = 0.;
			Real3D result = Real3D(0.,0.,0.);
			Real3D _myU = Real3D(0.,0.,0.);
			Real3D _ci = Real3D(0.,0.,0.);
			Int3D _Ni = latticeboltzmann->getMyNi();
			int _offset = latticeboltzmann->getHaloSkin();
			int _numVels = latticeboltzmann->getNumVels();
			
			for (int i = _offset; i < _Ni[0] - _offset; i++) {
				for (int j = _offset; j < _Ni[1] - _offset; j++) {
					for (int k = _offset; k < _Ni[2] - _offset; k++) {
						Real3D _jLoc = Real3D(0.,0.,0.);

						// calculation of density and momentum flux on the lattice site
						for (int l = 0; l < _numVels; l++) {
							_fi = latticeboltzmann->getLBFluid(Int3D(i,j,k),l);
							_ci = latticeboltzmann->getCi(l);
							_jLoc += _fi * _ci;
						}
                  // idea: maybe access LBMom directly and not compute them?
                  
						_myU += _jLoc;
					}
				}
			}
			_myU *= (latticeboltzmann->convTimeMDtoLB());
			_myU /= (latticeboltzmann->getA());
			
			boost::mpi::reduce(*getSystem()->comm, _myU, result, std::plus<Real3D>(), 0);

			if (getSystem()->comm->rank() == 0) {
				int _step = latticeboltzmann->getStepNum();
				printf ("statistics for step %d:\n", _step);
				
				printf("LB-fluid mom after streaming (LJ units) is %18.14f %18.14f %18.14f \n",
							 result.getItem(0), result.getItem(1), result.getItem(2));
			}
		}
		
		void LBOutputScreen::setLBTimerOld (real _lbTime_old) {lbTime_old = _lbTime_old;}
		real LBOutputScreen::getLBTimerOld() {return lbTime_old;}
		
		void LBOutputScreen::setLBTimerNew(real _lbTime_new) {lbTime_new = _lbTime_new;}
		real LBOutputScreen::getLBTimerNew() {return lbTime_new;}

		void LBOutputScreen::setOldStepNum(long int _oldStepNum) {oldStepNum = _oldStepNum;}
		long int LBOutputScreen::getOldStepNum() {return oldStepNum;}
		
    void LBOutputScreen::registerPython() {
      using namespace espressopp::python;

      class_<LBOutputScreen, bases< LBOutput > >
          ("analysis_LBOutput_Screen", init< shared_ptr< System >,
                                             shared_ptr< integrator::LatticeBoltzmann > >())

        .def("writeOutput", &LBOutputScreen::writeOutput)
      ;
    }
  }
}

