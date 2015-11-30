/*
  Copyright (C) 2012-2015
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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ConfigsParticleDecomp.hpp"
#include "bc/BC.hpp"
#include <boost/serialization/map.hpp>

using namespace std;
using namespace espressopp;

namespace espressopp {
  namespace analysis {

    using namespace iterator;

    int ConfigsParticleDecomp::getListSize() const{
      return configurations.size();
    }
   
    ConfigurationList ConfigsParticleDecomp::all() const{
      return configurations;
    }

    ConfigurationPtr ConfigsParticleDecomp::getConf(int position) const{
      int nconfigs = configurations.size();
      if (0 <= position and position < nconfigs) {
        return configurations[position];
      }
      else{
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        stringstream msg;
        msg << "Error. Velocities::get <out-of-range>" << endl;
        err.setException( msg.str() );
        return shared_ptr<Configuration>();
      }
    }

    void ConfigsParticleDecomp::pushConfig(ConfigurationPtr config){
      configurations.push_back(config);
    }
    
    void ConfigsParticleDecomp::gather() {
      System& system = getSystemRef();
      esutil::Error err(system.comm);
      
      int nprocs = system.comm->size();
      int myrank = system.comm->rank();
      
      int localN = system.storage->getNRealParticles();
      
      int curNumP = 0;
      boost::mpi::all_reduce(*system.comm, localN, curNumP, std::plus<int>());
      if(myrank==0){
        // check whether the number of particles is the same during the gathering
        if( curNumP != num_of_part ){
          stringstream msg;
          msg<<"   ConfigsParticleDecomp gathers the configurations of the same system\n"
                " with the same number of particles. If you need to store the systems\n"
                " with different number of particles you should use something else."
                " E.g `Configurations`";
          err.setException( msg.str() );
        }
      }

      ConfigurationPtr config = make_shared<Configuration> ();
      for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, Real3D > conf;
        if (rank_i == myrank) {
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            int id = cit->id();
            Real3D property = Real3D(0,0,0);
            if(key=="position")
              property = cit->position();
            else if(key=="velocity")
              property = cit->velocity();
            else if(key=="unfolded"){
              Real3D& pos = cit->position();
              Int3D& img = cit->image();
              Real3D Li = system.bc->getBoxL();
              for (int i = 0; i < 3; ++i) property[i] = pos[i] + img[i] * Li[i];
            }
            else{
              stringstream msg;
              msg<<"Error. Key "<<key<<" is unknown. Use position, unfolded or"
                      " velocity.";
              err.setException( msg.str() );
            }
            
            conf[id] = property;
          }
    	}

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        for (map<size_t,Real3D>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
          size_t id = itr->first;
          Real3D p = itr->second;
          if(idToCpu[id]==myrank) config->set(id, p[0], p[1], p[2]);
        }
      }

      pushConfig(config);
    }
    
    void ConfigsParticleDecomp::gatherFromFile(string filename) {
      System& system = getSystemRef();
      esutil::Error err(system.comm);

      int nprocs = system.comm->size();
      int myrank = system.comm->rank();

      int localN = system.storage->getNRealParticles();

      ConfigurationPtr config = make_shared<Configuration> ();
      map< size_t, Real3D > conf;

      if (myrank==0) {
          int id, type;
          real xpos, ypos, zpos;
          string line;
          ifstream file(filename.c_str());
          if (file.is_open()) {
        	// skip first 2 lines
          	getline(file, line);
        	getline(file, line);
        	int count = 0;
            while (getline(file, line)) {
              stringstream sl(line);
              sl >> id;
              sl >> type;
              sl >> xpos;
              sl >> ypos;
              sl >> zpos;
              // cout << id << ":" << x << "," << y << "," << z << endl;
              conf[id] = Real3D(xpos, ypos, zpos);
              count++;
            }
            file.close();
            cout << "read " << count << " particles from file " << filename << endl;
            if (count != num_of_part) {
                stringstream msg;
                msg << "Number of read particles does not match the number of particles of the system (which is " << num_of_part << ")";
                err.setException( msg.str() );
            }
          } else {
              stringstream msg;
              msg << "Unable to open file " << filename;
              err.setException( msg.str() );
          }
      }

      boost::mpi::broadcast(*system.comm, conf, 0);

      for (map<size_t,Real3D>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
        size_t id = itr->first;
        Real3D p = itr->second;
        if(idToCpu[id]==myrank) config->set(id, p[0], p[1], p[2]);
      }
      pushConfig(config);
    }

    // Python wrapping
    void ConfigsParticleDecomp::registerPython() {
      using namespace espressopp::python;

      class_<ConfigsParticleDecomp, boost::noncopyable >(
        "analysis_ConfigsParticleDecomp", no_init
        //init< shared_ptr< System > >()
      )
      .def_readonly("size", &ConfigsParticleDecomp::getListSize)
      
      .def("gather", &ConfigsParticleDecomp::gather)
      .def("gatherFromFile", &ConfigsParticleDecomp::gatherFromFile)
      .def("__getitem__", &ConfigsParticleDecomp::getConf)
      .def("all", &ConfigsParticleDecomp::all)
      .def("clear", &ConfigsParticleDecomp::clear)
      .def("compute", &ConfigsParticleDecomp::compute)
      ;
    }
  }
}
