/*
  Copyright (C) 2017
      Gregor Deichmann (TU Darmstadt, deichmann(at)cpc.tu-darmstadt.de) 
  
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

#include <fstream>
#include <iomanip>
#include "DumpXTC.hpp"
#include "storage/Storage.hpp"

#include "bc/BC.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

#include <boost/filesystem.hpp>

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;

namespace espressopp {
  namespace io {
    
    bool DumpXTC::open(const char *mode){

      if(mode[0] == 'a' && !boost::filesystem::exists(file_name)){
        fio = open_xtc(file_name.c_str(),"w"); //Opening with mode "a" on a non-existing file leads to error
      }
      else{
        fio = open_xtc(file_name.c_str(),mode); 
      }

      return true;
    }

    void DumpXTC::close(){

      close_xtc(fio);
      
      return;
    }
      
    void DumpXTC::dump(){

      shared_ptr<System> system = getSystem();
      ConfigurationsExt conf( system );
      conf.setUnfolded(unfolded);
      conf.gather();
      
      if( system->comm->rank()==0 ){
        ConfigurationExtPtr conf_real = conf.back();
        
        int num_of_particles = conf_real->getSize();
       
        if(this->open("a")){ 
        
          ConfigurationExtIterator cei = conf_real-> getIterator();
          //Real3D *box = new Real3D [dim];
          matrix box;
          rvec *coord = new rvec [num_of_particles];
          RealND props;
          props.setDimension( cei.currentProperties().getDimension() );

          //HACK: Only valid for orthorhombic BC. Will there be anything else in ESPP?
          Real3D bl=system->bc->getBoxL();

          for(int i=0;i<dim;i++){
            box[i][0] = 0.;
            box[i][1] = 0.;
            box[i][2] = 0.;

            box[i][i] = bl[i]; 
          }

          for(int i=0;i<num_of_particles;i++){
            
            props=cei.nextProperties();
            coord[i][0] = props[0]*length_factor; //We only write coordinates to .xtc 
            coord[i][1] = props[1]*length_factor;
            coord[i][2] = props[2]*length_factor;

          }

          int step = integrator->getStep();
          float time = integrator->getTimeStep()*step;

          write_xtc(fio, num_of_particles, step, time, box, coord, xtcprec);

          delete [] coord;
          
          this->close();
        }
        else cout << "Unable to open file: "<< file_name <<endl;

      
      }

      return;
    }
      
    // Python wrapping
    void DumpXTC::registerPython() {

      using namespace espressopp::python;

      class_<DumpXTC, bases<ParticleAccess>, boost::noncopyable >
      ("io_DumpXTC", init< shared_ptr< System >, 
                           shared_ptr< integrator::MDIntegrator >, 
                           std::string, 
                           bool,
                           real,
                           bool>())
        .add_property("filename", &DumpXTC::getFilename, 
                                  &DumpXTC::setFilename)
        .add_property("unfolded", &DumpXTC::getUnfolded, 
                                  &DumpXTC::setUnfolded)
        .add_property("append",   &DumpXTC::getAppend, 
                                  &DumpXTC::setAppend)
        .def("dump", &DumpXTC::dump)
      ;
    }
  }
}
