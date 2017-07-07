/*
  Copyright (C) 2017
      Gregor Deichmann (TU Darmstadt, deichmann(at)cpc.tu-darmstadt.de) 
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

#include <fstream>
#include <iomanip>
#include "DumpXTC.hpp"
#include "storage/Storage.hpp"

#include "bc/BC.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

extern "C"{
#include "libxdrf.c" //c-routines for xdr io and compression
}

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;

namespace espressopp {
  namespace io {
    
    bool DumpXTC::open(const char *mode){
      
      fp = fopen(file_name.c_str(),mode);
      
      if(!fp){
        return false;
      }
      
      xdrstdio_create(xdr,fp,xdrmode);
      
      return true;
    }

    void DumpXTC::close(){
      
      if(!fp){
        return;
      }
      
      if(xdr){
        xdr_destroy(xdr);
      }
      
      if(fp){
        fclose(fp);
        fp=NULL;
      }   
      
      return;

    }

    bool DumpXTC::write(int natoms,  
                              int step,
                              float time,
                              Real3D *box,
                              float *x,
                              float prec){

        bool bOk = true;

        int num=XTC_MAGIC,ival; //gromacs magic number (should be 1995)
        float fval;
        int i,j;

        if( xdr_int(xdr,&num) == 0){
            bOk = false;
            return 1;
        }

        ival=natoms;
        if( xdr_int(xdr,&ival) == 0){
            bOk = false;
            return 1;
        }

        ival=step;
        if( xdr_int(xdr,&ival) == 0){
            bOk = false;
            return 1;
        };

        fval=time;
        if( xdr_float(xdr,&fval) == 0){
            bOk = false;
            return 1;
        };

        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
                fval = box[i][j];

                if( xdr_float(xdr,&fval) == 0){
                    bOk = false;
                    break;
                }
            }
            if(!bOk){
                break;
            }
        }

        if(!bOk){
            return 1;
        }

        const int x_size = natoms*dim;
        
        if( xdr3dfcoord(xdr, x, &natoms, &prec) == 0 ){
            bOk=false;
        }

        return bOk;
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
          //RealND::iterator ii;
          //short ind;
          Real3D *box = new Real3D [dim];
          float *coord = new float [dim*num_of_particles];
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
            coord[i*dim] = props[0]*length_factor; //We only write coordinates to .xtc 
            coord[i*dim+1] = props[1]*length_factor;
            coord[i*dim+2] = props[2]*length_factor;

          }

          int step = integrator->getStep();
          float time = integrator->getTimeStep()*step;

          this->write(num_of_particles,step,time,box,coord,xtcprec);

          delete [] box;
          delete [] coord;
          
          this->close();
        }
        else cout << "Unable to open file: "<< file_name <<endl;

      
      }
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
                           std::string,
                           bool>())
        .add_property("filename", &DumpXTC::getFilename, 
                                  &DumpXTC::setFilename)
        .add_property("unfolded", &DumpXTC::getUnfolded, 
                                  &DumpXTC::setUnfolded)
        .add_property("append", &DumpXTC::getAppend, 
                                  &DumpXTC::setAppend)
        .add_property("length_factor", &DumpXTC::getLengthFactor, 
                                       &DumpXTC::setLengthFactor)
        .add_property("length_unit", &DumpXTC::getLengthUnit, 
                                     &DumpXTC::setLengthUnit)
        .def("dump", &DumpXTC::dump)
      ;
    }
  }
}
