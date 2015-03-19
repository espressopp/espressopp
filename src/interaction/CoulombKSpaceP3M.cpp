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
//#include <boost/signals2.hpp>
#include "CoulombKSpaceP3M.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class CellListAllParticlesInteractionTemplate <CoulombKSpaceP3M>
    CellListCoulombKSpaceP3M;

    CoulombKSpaceP3M::
    CoulombKSpaceP3M(shared_ptr< System > _system,
                     real _coulomb_prefactor,
                     real _alpha,
                     Int3D _M,
                     int _P,
                     real _rcut,
                     int _interpolation
              ): system(_system), C_pref(_coulomb_prefactor), alpha(_alpha),
                    M(_M), P(_P), rc(_rcut), interpolation(_interpolation){
      
      // predefined assigned function coefficients
      af_coef[1][0][0] = 1.0;

      af_coef[2][0][0] =   0.5;
      af_coef[2][0][1] =  -1.0;
      af_coef[2][1][0] =   0.5;
      af_coef[2][1][1] =  -1.0;

      af_coef[3][0][0] =   1./8.;
      af_coef[3][0][1] =  -4./8.;
      af_coef[3][0][2] =   4./8.;
      af_coef[3][1][0] =   6./8.;
      af_coef[3][1][1] =   0.;
      af_coef[3][1][2] =   8./8.;
      af_coef[3][2][0] =   1./8.;
      af_coef[3][2][1] =   4./8.;
      af_coef[3][2][2] =   4./8.;


      af_coef[4][0][0] =   1./48.;
      af_coef[4][0][1] = - 6./48.;
      af_coef[4][0][2] =  12./48.;
      af_coef[4][0][3] = - 8./48.;
      af_coef[4][1][0] =  23./48.;
      af_coef[4][1][1] = -30./48.;
      af_coef[4][1][2] = -12./48.;
      af_coef[4][1][3] =  24./48.;
      af_coef[4][2][0] =  23./48.;
      af_coef[4][2][1] =  30./48.;
      af_coef[4][2][2] = -12./48.;
      af_coef[4][2][3] = -24./48.;
      af_coef[4][0][0] =   1./48.;
      af_coef[4][0][1] = - 6./48.;
      af_coef[4][0][2] =  12./48.;
      af_coef[4][0][3] = - 8./48.;

      af_coef[5][0][0] =    1./384.;
      af_coef[5][0][1] =  - 8./384.;
      af_coef[5][0][2] =   24./384.;
      af_coef[5][0][3] =  -32./384.;
      af_coef[5][0][4] =   16./384.;
      af_coef[5][1][0] =   76./384.;
      af_coef[5][1][1] = -176./384.;
      af_coef[5][1][2] =   96./384.;
      af_coef[5][1][3] =   64./384.;
      af_coef[5][1][4] =  -64./384.;
      af_coef[5][2][0] =  230./384.;
      af_coef[5][2][1] =        0.0;
      af_coef[5][2][2] = -240./384.;
      af_coef[5][2][3] =        0.0;
      af_coef[5][2][4] =   96./384.;
      af_coef[5][3][0] =   76./384.;
      af_coef[5][3][1] =  176./384.;
      af_coef[5][3][2] =   96./384.;
      af_coef[5][3][3] =  -64./384.;
      af_coef[5][3][4] =  -64./384.;
      af_coef[5][4][0] =    1./384.;
      af_coef[5][4][1] =    8./384.;
      af_coef[5][4][2] =   24./384.;
      af_coef[5][4][3] =   32./384.;
      af_coef[5][4][4] =   16./384.;

      af_coef[6][0][0] =     1./3840.;
      af_coef[6][0][1] =   -10./3840.;
      af_coef[6][0][2] =    40./3840.;
      af_coef[6][0][3] =   -80./3840.;
      af_coef[6][0][4] =    80./3840.;
      af_coef[6][0][5] =    32./3840.;
      af_coef[6][1][0] =   237./3840.;
      af_coef[6][1][1] =  -750./3840.;
      af_coef[6][1][2] =   840./3840.;
      af_coef[6][1][3] =  -240./3840.;
      af_coef[6][1][4] =  -240./3840.;
      af_coef[6][1][5] =   160./3840.;
      af_coef[6][2][0] =  1682./3840.;
      af_coef[6][2][1] = -1540./3840.;
      af_coef[6][2][2] =  -880./3840.;
      af_coef[6][2][3] =  1120./3840.;
      af_coef[6][2][4] =   160./3840.;
      af_coef[6][2][5] =  -320./3840.;
      af_coef[6][3][0] =  1682./3840.;
      af_coef[6][3][1] =  1540./3840.;
      af_coef[6][3][2] =  -880./3840.;
      af_coef[6][3][3] = -1120./3840.;
      af_coef[6][3][4] =   160./3840.;
      af_coef[6][3][5] =   320./3840.;
      af_coef[6][4][0] =   237./3840.;
      af_coef[6][4][1] =   750./3840.;
      af_coef[6][4][2] =   840./3840.;
      af_coef[6][4][3] =   240./3840.;
      af_coef[6][4][4] =  -240./3840.;
      af_coef[6][4][5] =  -160./3840.;
      af_coef[6][5][0] =     1./3840.;
      af_coef[6][5][1] =    10./3840.;
      af_coef[6][5][2] =    40./3840.;
      af_coef[6][5][3] =    80./3840.;
      af_coef[6][5][4] =    80./3840.;
      af_coef[6][5][5] =    32./3840.;

      af_coef[7][0][0] =      1./46080.;
      af_coef[7][0][1] =    -12./46080.;
      af_coef[7][0][2] =     60./46080.;
      af_coef[7][0][3] =   -160./46080.;
      af_coef[7][0][4] =    240./46080.;
      af_coef[7][0][5] =   -192./46080.;
      af_coef[7][0][6] =     64./46080.;
      af_coef[7][1][0] =    722./46080.;
      af_coef[7][1][1] =  -2832./46080.;
      af_coef[7][1][2] =   4440./46080.;
      af_coef[7][1][3] =  -3200./46080.;
      af_coef[7][1][4] =    480./46080.;
      af_coef[7][1][5] =    768./46080.;
      af_coef[7][1][6] =   -384./46080.;
      af_coef[7][2][0] =  10543./46080.;
      af_coef[7][2][1] = -17340./46080.;
      af_coef[7][2][2] =   4740./46080.;
      af_coef[7][2][3] =   6880./46080.;
      af_coef[7][2][4] =  -4080./46080.;
      af_coef[7][2][5] =   -960./46080.;
      af_coef[7][2][6] =    960./46080.;
      af_coef[7][3][0] =  23548./46080.;
      af_coef[7][3][1] =            0.0;
      af_coef[7][3][2] = -18480./46080.;
      af_coef[7][3][3] =            0.0;
      af_coef[7][3][4] =   6720./46080.;
      af_coef[7][3][5] =            0.0;
      af_coef[7][3][6] =  -1280./46080.;
      af_coef[7][4][0] =  10543./46080.;
      af_coef[7][4][1] =  17340./46080.;
      af_coef[7][4][2] =   4740./46080.;
      af_coef[7][4][3] =  -6880./46080.;
      af_coef[7][4][4] =  -4080./46080.;
      af_coef[7][4][5] =    960./46080.;
      af_coef[7][4][6] =    960./46080.;
      af_coef[7][5][0] =    722./46080.;
      af_coef[7][5][1] =   2832./46080.;
      af_coef[7][5][2] =   4440./46080.;
      af_coef[7][5][3] =   3200./46080.;
      af_coef[7][5][4] =    480./46080.;
      af_coef[7][5][5] =   -768./46080.;
      af_coef[7][5][6] =   -384./46080.;
      af_coef[7][6][0] =      1./46080.;
      af_coef[7][6][1] =     12./46080.;
      af_coef[7][6][2] =     60./46080.;
      af_coef[7][6][3] =    160./46080.;
      af_coef[7][6][4] =    240./46080.;
      af_coef[7][6][5] =    192./46080.;
      af_coef[7][6][6] =     64./46080.;
      
      getParticleNumber();
      preset();
      
      // first initialize
      initialize();
        
      // This function calculates the square of all particle charges. It should be called ones,
      // if the total number of particles doesn't change.
      count_charges(system->storage->getRealCells()); 

      /* make a connection to boundary conditions to invoke recalculation of KVec if box
         dimensions change
      */
      connectionRecalcKVec = system->bc->onBoxDimensionsChanged.
              connect(boost::bind(&CoulombKSpaceP3M::preset, this));
      // make a connection to storage in order to get number of particles
      connectionGetParticleNumber = system->storage->onParticlesChanged.
              connect(boost::bind(&CoulombKSpaceP3M::getParticleNumber, this));
      
      // sign a signal in order not to recalculate common part twice
      //recalcCommonPart = system.
    }
    
    CoulombKSpaceP3M::~CoulombKSpaceP3M(){
      /*
      delete [] sum;
      delete [] totsum;
      sum = NULL;
      totsum = NULL;
      */
      clean_fftw();
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void CoulombKSpaceP3M::registerPython() {
      using namespace espressopp::python;

      class_< CoulombKSpaceP3M, bases< Potential > >
      ("interaction_CoulombKSpaceP3M", 
              init< shared_ptr<System>, real, real, Int3D, int, real, int >() )
    	.add_property("prefactor", &CoulombKSpaceP3M::getPrefactor, 
                                   &CoulombKSpaceP3M::setPrefactor);
    	//.add_property("alpha", &CoulombKSpaceP3M::getAlpha, &CoulombKSpaceP3M::setAlpha)
    	//.add_property("kmax", &CoulombKSpaceP3M::getKMax, &CoulombKSpaceP3M::setKMax)
      //;

      class_< CellListCoulombKSpaceP3M, bases< Interaction > >
        ("interaction_CellListCoulombKSpaceP3M",
              init< shared_ptr< storage::Storage >,
                    shared_ptr< CoulombKSpaceP3M > >())
        .def("getPotential", &CellListCoulombKSpaceP3M::getPotential)
	  ;

    }
    
  }
}
