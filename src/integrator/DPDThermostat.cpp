/*
  
  Copyright (C) 2017
      Gregor Deichmann (TU Darmstadt)
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
#include "DPDThermostat.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;

    DPDThermostat::DPDThermostat(
    		shared_ptr<System> system,
    		shared_ptr<VerletList> _verletList)
    :Extension(system), verletList(_verletList) {

      type = Extension::Thermostat;

      gamma  = 0.0;
      temperature = 0.0;

      current_cutoff = verletList->getVerletCutoff() - system->getSkin();
      current_cutoff_sqr = current_cutoff*current_cutoff;
      
      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "DPD constructed");
    }

    void DPDThermostat::setGamma(real _gamma)
    {
      gamma = _gamma;
    }

    real DPDThermostat::getGamma()
    {
      return gamma;
    }

    void DPDThermostat::setTGamma(real _tgamma) {
      tgamma = _tgamma;
    }

    real DPDThermostat::getTGamma() {
      return tgamma;
    }

    void DPDThermostat::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real DPDThermostat::getTemperature()
    {
      return temperature;
    }

    DPDThermostat::~DPDThermostat() {
        disconnect();
    }


    void DPDThermostat::disconnect() {

        _initialize.disconnect();
        _heatUp.disconnect();
        _coolDown.disconnect();
        _thermalize.disconnect();
    }

    void DPDThermostat::connect() {

        // connect to initialization inside run()
        _initialize = integrator->runInit.connect(
                boost::bind(&DPDThermostat::initialize, this));

        _heatUp = integrator->recalc1.connect(
                boost::bind(&DPDThermostat::heatUp, this));

        _coolDown = integrator->recalc2.connect(
                boost::bind(&DPDThermostat::coolDown, this));

        _thermalize = integrator->aftInitF.connect(
                boost::bind(&DPDThermostat::thermalize, this));
    }


    void DPDThermostat::thermalize() {
        LOG4ESPP_DEBUG(theLogger, "thermalize DPD");

        System& system = getSystemRef();
        system.storage->updateGhostsV();

        // loop over VL pairs
        for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
            Particle &p1 = *it->first;
            Particle &p2 = *it->second;
      
            if(gamma > 0.0)
              frictionThermoDPD(p1, p2);
            if(tgamma > 0.0)
              frictionThermoTDPD(p1, p2);
        }
    }


    void DPDThermostat::frictionThermoDPD(Particle& p1, Particle& p2) {
      //Implements the standard DPD thermostat 
      Real3D r = p1.position() - p2.position();
      real dist2 = r.sqr();
      
      if(dist2 < current_cutoff_sqr){
        
        real dist = sqrt(dist2);
        real omega = 1-dist/current_cutoff;
        real omega2 = omega*omega;

        r /= dist;

        real veldiff = (p1.velocity() - p2.velocity()) * r;
        real friction = pref1 * omega2 * veldiff;
        real r0 = ((*rng)() - 0.5);
        real noise = pref2 * omega * r0;//(*rng)() - 0.5);

        Real3D f = (noise - friction) * r;
        p1.force() += f;
        p2.force() -= f;
      }
    }

    void DPDThermostat::frictionThermoTDPD(Particle& p1, Particle& p2) {
      //Implements a transverse DPD thermostat with the canonical functional form of omega
      Real3D r = p1.position() - p2.position();
	  real dist2 = r.sqr();
      
      if(dist2 < current_cutoff_sqr){
        
        real dist = sqrt(dist2);
        real omega = 1-dist/current_cutoff;
        real omega2 = omega*omega;

        r /= dist;
		
        Real3D noisevec(0.0);
        noisevec[0] = (*rng)() - 0.5;
        noisevec[1] = (*rng)() - 0.5;
        noisevec[2] = (*rng)() - 0.5;
        
        Real3D veldiff = p1.velocity() - p2.velocity();

        Real3D f_damp,f_rand;
        
        //Calculate matrix product of projector and veldiff vector:
        //P dv = (I - r r_T) dv 
        f_damp[0] = (1.0 - r[0]*r[0])*veldiff[0] - r[0]*r[1]*veldiff[1] - r[0]*r[2]*veldiff[2];
        f_damp[1] = (1.0 - r[1]*r[1])*veldiff[1] - r[1]*r[0]*veldiff[0] - r[1]*r[2]*veldiff[2];
        f_damp[2] = (1.0 - r[2]*r[2])*veldiff[2] - r[2]*r[0]*veldiff[0] - r[2]*r[1]*veldiff[1];
       
        //Same with random vector
        f_rand[0] = (1.0 - r[0]*r[0])*noisevec[0] - r[0]*r[1]*noisevec[1] - r[0]*r[2]*noisevec[2];
        f_rand[1] = (1.0 - r[1]*r[1])*noisevec[1] - r[1]*r[0]*noisevec[0] - r[1]*r[2]*noisevec[2];
        f_rand[2] = (1.0 - r[2]*r[2])*noisevec[2] - r[2]*r[0]*noisevec[0] - r[2]*r[1]*noisevec[1];
        
        f_damp *= pref3 * omega2;
        f_rand *= pref4 * omega;

		p1.force() += f_rand - f_damp;
		p2.force() -= f_rand - f_damp;
      }	
	}

    void DPDThermostat::initialize() {
    	// calculate the prefactors
      System& system = getSystemRef();
      current_cutoff = verletList->getVerletCutoff() - system.getSkin();
      current_cutoff_sqr = current_cutoff*current_cutoff;
      
      real timestep = integrator->getTimeStep();

      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
              ", gamma = " << gamma <<
              ", tgamma = " << tgamma <<
              ", temperature = " << temperature);

      pref1 = gamma;
      pref2 = sqrt(24.0 * temperature * gamma / timestep);
      pref3 = tgamma;
      pref4 = sqrt(24.0 * temperature * tgamma / timestep);
    }

    /** very nasty: if we recalculate force when leaving/reentering the integrator,
	a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
	numbers are drawn twice, resulting in a different variance of the random force.
	This is corrected by additional heat when restarting the integrator here.
	Currently only works for the Langevin thermostat, although probably also others
	are affected.
    */

    void DPDThermostat::heatUp() {
    	LOG4ESPP_INFO(theLogger, "heatUp");

    	
        pref2buffer = pref2;
    	pref2       *= sqrt(3.0);
    	pref4buffer = pref4;
    	pref4       *= sqrt(3.0);
        
    }

    /** Opposite to heatUp */

    void DPDThermostat::coolDown(){
        LOG4ESPP_INFO(theLogger, "coolDown");

        
        pref2 = pref2buffer;
        pref4 = pref4buffer;
        
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void DPDThermostat::registerPython() {
      using namespace espressopp::python;
      class_<DPDThermostat, shared_ptr<DPDThermostat>, bases<Extension> >
        ("integrator_DPDThermostat", init<shared_ptr<System>, shared_ptr<VerletList> >())
        .def("connect", &DPDThermostat::connect)
        .def("disconnect", &DPDThermostat::disconnect)
        .add_property("gamma", &DPDThermostat::getGamma, &DPDThermostat::setGamma)
        .add_property("tgamma", &DPDThermostat::getTGamma, &DPDThermostat::setTGamma)
        .add_property("temperature", &DPDThermostat::getTemperature, &DPDThermostat::setTemperature)
        ;
    }
  }
}

