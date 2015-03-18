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
#include "BerendsenBarostat.hpp"
#include "System.hpp"
#include "esutil/Error.hpp"

namespace espressopp {
  
  using namespace analysis;
  using namespace esutil;
  using namespace std;
  
  namespace integrator {

    LOG4ESPP_LOGGER(BerendsenBarostat::theLogger, "BerendsenBarostat");

    BerendsenBarostat::BerendsenBarostat(shared_ptr<System> system): Extension(system){
      tau  = 1.0;
      P0 = 1.0;
      
      fixed = Int3D(1,1,1);
      
      exponent = 1./3.;
      
      type = Extension::Barostat;

      LOG4ESPP_INFO(theLogger, "BerendsenBarostat constructed");
    }

    BerendsenBarostat::~BerendsenBarostat(){
      LOG4ESPP_INFO(theLogger, "~BerendsenBarostat");
      disconnect();
    }
    
    void BerendsenBarostat::disconnect(){
      _runInit.disconnect();
      _aftIntV.disconnect();
    }

    void BerendsenBarostat::connect(){
      // connection to initialisation
      _runInit = integrator->runInit.connect( boost::bind(&BerendsenBarostat::initialize, this));

      // connection to the signal at the end of the run
      _aftIntV = integrator->aftIntV.connect( boost::bind(&BerendsenBarostat::barostat, this));
    }

    // set and get time constant for Berendsen barostat
    void BerendsenBarostat::setTau(real _tau) {
      tau = _tau;
    }
    real BerendsenBarostat::getTau() {
      return tau;
    }
    // set and get external pressure
    void BerendsenBarostat::setPressure(real _P0) {
      P0 = _P0;
    }
    real BerendsenBarostat::getPressure() {
      return P0;
    }
    void BerendsenBarostat::setFixed(Int3D _fix) {
      fixed = _fix;
      int e=0;
      for(int i=0;i<3;i++)
        e+=fixed[i];
      exponent = 1.0/(real)e;
    }
    Int3D BerendsenBarostat::getFixed() {
      return fixed;
    }

    void BerendsenBarostat::barostat(){
      LOG4ESPP_DEBUG(theLogger, "equilibrating the pressure");

      System& system = getSystemRef();
      
      Pressure Pcurrent(getSystem());
      
      real P = Pcurrent.compute();  // calculating the current pressure in system
      
      real mu3 = 1 + pref * (P - P0);
      
      Error err(system.comm);
      if(mu3<0.0){
        stringstream msg;
        msg << "Scaling coefficient is <0 (Berendsen barostat). mu^3="<<mu3;
        err.setException( msg.str() );
        err.checkException();
      }
      
      real mu = pow(  mu3, exponent );  // calculating the current scaling parameter
      
      Real3D mu3D( 1.0 );
      for(int i=0;i<3;i++)
        if( fixed[i]>0 )
          mu3D[i]=mu;

      system.scaleVolume( mu3D, true );
    }

     // calculate the prefactors
    void BerendsenBarostat::initialize(){
      LOG4ESPP_INFO(theLogger, "init, tau = " << tau << 
                               ", external pressure = " << P0);
      real dt = integrator->getTimeStep();
      pref = dt / tau;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void BerendsenBarostat::registerPython() {

      using namespace espressopp::python;

      class_<BerendsenBarostat, shared_ptr<BerendsenBarostat>, bases<Extension> >

        ("integrator_BerendsenBarostat", init< shared_ptr<System> >())

        .add_property("tau",
              &BerendsenBarostat::getTau,
              &BerendsenBarostat::setTau)
        .add_property("pressure",
              &BerendsenBarostat::getPressure, 
              &BerendsenBarostat::setPressure)
        .add_property("fixed",
              &BerendsenBarostat::getFixed, 
              &BerendsenBarostat::setFixed)
      
        .def("connect", &BerendsenBarostat::connect)
        .def("disconnect", &BerendsenBarostat::disconnect)
      ;
    }

  }
}

