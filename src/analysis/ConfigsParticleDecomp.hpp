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

namespace espressopp {
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
      
       /*
       * Constructor, allow for unlimited snapshots. It defines how many particles
       * correspond to different cpu without breaking chains. So the monomers of
       * one chain correspond to one CPU only.
        * 
        * !! currently only works for particles numbered like 0, 1, 2,... !!
        * !! with each chain consisting particles with subsequent ids     !!
       */
      ConfigsParticleDecomp(shared_ptr<System> system, int _chainlength): SystemAccess (system){
        // by default key = "position", it will store the particle positions
        // (option: "velocity" or "unfolded")
        esutil::Error err(system->comm);
        
        key = "position";
        chainlength = _chainlength;
        
        int localN = system -> storage -> getNRealParticles();
        boost::mpi::all_reduce(*system->comm, localN, num_of_part, std::plus<int>());
        
        int n_nodes = system -> comm -> size();
        int this_node = system -> comm -> rank();
              
        //for monodisperse chains
        int num_chains = num_of_part / chainlength;
        int local_num_chains = (int) ceil( (double)num_chains / n_nodes );
        int local_num_part = local_num_chains * chainlength;
        
        //in case the chainlength does not match the total number of particles
        if(num_of_part % chainlength != 0){
            cout << "chainlength does not match total number of particles\n"
                    << "chainlength: " << chainlength
                    << "\n num_of_part " << num_of_part << "\n\n";             
        }        
        
        //assignment particles to cpus (= filling of map idToCpu)
        //CPU0 will use particles 0, 1, 2, ... local_num_particles-1.
        //CPU1 will use particles local_num_particles, local_num_particles+1,...        
        int nodeNum = -1;
        for(long unsigned int id = 0; id < num_of_part ;id++){
            if(id % local_num_part == 0){ ++nodeNum;}
            idToCpu[id] = nodeNum;
            //intId++   //...if loop was with iterator
            //if(intId %)local_num_part == 0) nodeNum++;            
        }
        //output if the assignment failed
        if (nodeNum >= n_nodes) {
            if(this_node == 0){
                cout << "assignment went wrong. Particles were assigned to proc "
                        << nodeNum << "\n";
                cout << "highest process number should be " << n_nodes - 1 <<"\n";
                cout << "check if total number of particles matches with chainlength\n";
            }
        }
        
                
        //connecting particle ids with their chain ids (= filling of map idToCid)
        //chain 0 consists of particles 0, 1, 2, ... chainlength -1
        //chain1 consists of particles chainlength, chainlength+1,...        
        int cid = -1;
        for(long unsigned int id = 0; id < num_of_part ;id++){
            if(id % chainlength == 0){ ++cid;}
            idToCid[id] = cid;           
        }
        //output if the assignment failed
        if (cid >= num_chains) {
            if(this_node == 0){
                cout << "assignment of chain ids went wrong."
                     << "Particles were assigned to chain "
                        << cid << "\n";
                cout << "highest chain id should be " << num_chains - 1 <<"\n";
                cout << "check if total number of particles matches with chainlength\n";
            }
        }
      
        
        //assignment chains to cpus (= filling of map cidToCpu)
        //CPU0 will use particles 0, 1, 2, ... local_num_particles-1.
        //CPU1 will use particles local_num_particles, local_num_particles+1,...        
        nodeNum = -1;
        //using k instead of cid to avoid conflicts with global cid
        for(long unsigned int k = 0; k < num_chains ;k++){
            if(k % local_num_chains == 0){ ++nodeNum;}
            cidToCpu[k] = nodeNum;
            //intId++   //...if loop was with iterator
            //if(intId %)local_num_part == 0) nodeNum++;            
        }
        //output if the assignment failed
        if (nodeNum >= n_nodes) {
            if(this_node == 0){
                cout << "assignment went wrong. Chains were assigned to proc "
                        << nodeNum << "\n";
                cout << "highest process number should be " << n_nodes - 1 <<"\n";
                cout << "check if total number of particles matches with chainlength\n";
            }
        } 
        
        try{
          if(num_chains < n_nodes){
            stringstream msg;
            msg<<"Warning. Number of chains less then the number of nodes.\n";
            msg<<"It might be a problem. NChains="<<num_chains<<" NNodes="<<n_nodes;
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
      int chainlength; //for calculations with chains (instead of monomers)
      map< size_t, int > idToCpu; // binds cpu and particle id
      map< size_t, int > cidToCpu; // binds cpu and chain id
      map< size_t, int > idToCid; // binds cpu and particle id
      
      string key; // it can be "position", "velocity" or "unfolded"
      
    private:

      void pushConfig(ConfigurationPtr config);
 
      // the list of snapshots
      ConfigurationList configurations;
    };
  }
}

#endif
