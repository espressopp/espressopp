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

// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGSPARTICLEDECOMP_HPP
#define _ANALYSIS_CONFIGSPARTICLEDECOMP_HPP

#include "python.hpp"
#include "mpi.h"
#include "types.hpp"
#include "SystemAccess.hpp"
#include "Configuration.hpp"

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"

#include <string>

using namespace std;

namespace espresso {
  namespace analysis {

    using namespace iterator;
    /*
     * Class that stores particle !!properties (velocities at the moment)!! for later
     * analysis. It uses object Configuration to store data.
     * 
     * Here the concept of particle decomposition is used, i.e. each processor stores
     * relevant number of particles. It's useless to get the data on python level from
     * here. Therefore it is abstract class. A derived class should realize the function
     * `compute`.
     * 
     * Important: Mainly it was created in order to observe the system in time.
     * !!At the moment the number of particles should be the same for different snapshots.!!
     * Otherwise it will throw a runtime error exception
    */

    typedef vector<ConfigurationPtr> ConfigurationList;

    class ConfigsParticleDecomp : public SystemAccess {

    public:
      /*
       * Constructor, allow for unlimited snapshots. It defines how many particles
       * correspond to different cpu.
       */
      ConfigsParticleDecomp(shared_ptr<System> system): SystemAccess (system){
        // by default key = "position", it will store the particle positions
        // (option: "velocity" or "unfolded")
        esutil::Error err(system->comm);
        
        key = "position";
        
        int localN = system -> storage -> getNRealParticles();
        boost::mpi::all_reduce(*system->comm, localN, num_of_part, std::plus<int>());
        
        int n_nodes = system -> comm -> size();
        int this_node = system -> comm -> rank();
        
        int local_num_of_part = num_of_part / n_nodes + 1;

        vector<int> tot_idList;
        for(int rank_i = 0; rank_i<n_nodes;rank_i++){
          
          int numLocPart = 0;
          if(rank_i==this_node){
            numLocPart = system -> storage -> getNRealParticles();
          }
          boost::mpi::broadcast(*system->comm, numLocPart, rank_i);
          
          int* idList = new int[numLocPart];
          
          if(rank_i==this_node){
            
            int count = 0;
            CellList realCells = system -> storage -> getRealCells();
            for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
              int id = cit->id();

              idList[count] = id;
              count++;
            }
          }
          
          boost::mpi::broadcast(*system->comm, idList, numLocPart, rank_i);
          
          for(int i=0; i<numLocPart;i++){
            tot_idList.push_back( idList[i] );
          }
          
          delete [] idList;
          idList = NULL;
        }

        int nodeNum = 0;
        int count = 0;
        for (vector<int>::iterator it = tot_idList.begin(); it!=tot_idList.end(); ++it) {
          idToCpu[*it] = nodeNum;
          count ++;
          if(count>=local_num_of_part){
            count = 0;
            nodeNum++;
          }
        }        
        
        try{
          if(num_of_part <= n_nodes){
            stringstream msg;
            msg<<"Warning. Number of particles less then the number of nodes.\n";
            msg<<"It might be a problem. NPart="<<num_of_part<<" NNodes="<<n_nodes;
            err.setException( msg.str() );
            err.checkException();
          }
        }
        catch(std::exception const& e){
          if(this_node==0)
            cout << "Exception: " << e.what() << "\n";
        }
      }
      ~ConfigsParticleDecomp() {}

      // get number of available snapshots. Returns the size of Configurationlist
      int getListSize() const;

      // Take a snapshot of property (all current particle velocities at the moment)
      void gather();

      // Read in a snapshot from a xyz-file
      void gatherFromFile(string filename);

      // Get a configuration from ConfigurationList
      ConfigurationPtr getConf(int position) const;

      // it returns all the configurations
      ConfigurationList all() const;

      // it erases all the configurations from ConfigurationList
      void clear(){
        configurations.clear();
      }
      
      virtual python::list compute() const = 0;

      static void registerPython();
    
    protected:

      static LOG4ESPP_DECL_LOGGER(logger);

      // all cpus handle defined number of particles
      int num_of_part;
      
      map< size_t, int > idToCpu; // binds cpu and particle id
      
      string key; // it can be "position", "velocity" or "unfolded"
      
    private:

      void pushConfig(ConfigurationPtr config);
 
      // the list of snapshots
      ConfigurationList configurations;
    };
  }
}

#endif
