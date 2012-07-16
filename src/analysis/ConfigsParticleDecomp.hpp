// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGSPARTICLEDECOMP_HPP
#define _ANALYSIS_CONFIGSPARTICLEDECOMP_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "Configuration.hpp"

#include "storage/Storage.hpp"
#include "python.hpp"

#include <string>

using namespace std;

namespace espresso {
  namespace analysis {

    /*
     * Class that stores particle !!properties (velocities at the moment)!! for later
     * analysis. It uses object Configuration to store data.
     * 
     * Here the concept of particle decomposition is used, i.e. each processor stores
     * relevant number of particles. It's useless to get the data on python level from
     * here. Therefore it is abstract class. A derived class should realize the function
     * `calculate`.
     * 
     * Important: Mainly it was created in order to observe the system in time.
     * At the moment the number of particles should be the same for different snapshots.
    */

    typedef vector<ConfigurationPtr> ConfigurationList;

    class ConfigsParticleDecomp : public SystemAccess {

    public:
      // Constructor, allow for unlimited snapshots. It defines how many particles
      // correspond to different cpu.
      // TODO !!Warning. Now it works only for sequential id 0 - num_of_part
      ConfigsParticleDecomp(shared_ptr<System> system): SystemAccess (system){
        // by default key = "position", it will store the particle positions
        key = "position";
        
        int localN = system -> storage -> getNRealParticles();
        boost::mpi::all_reduce(*system->comm, localN, num_of_part, std::plus<int>());
        
        int n_nodes = system -> comm -> size();
        int this_node = system -> comm -> rank();
        
        // TODO it could be a problem if   n_nodes > kVectorLength !!!
        local_num_of_part = num_of_part / n_nodes + 1;

        min_id = this_node * local_num_of_part;
        max_id = min_id + local_num_of_part;
        if( max_id>num_of_part ) max_id = num_of_part;
        
        if(num_of_part <= n_nodes){
          cout<<"Warning. Number of particles less then the number of nodes."<<endl;
          cout<<"nparts="<<num_of_part<<"  nnodes="<<n_nodes<<endl;
        }
      }
      ~ConfigsParticleDecomp() {}

      // get number of available snapshots. Returns the size of Configurationlist
      int getListSize() const;

      // Take a snapshot of property (all current particle velocities at the moment)
      void gather();

      // Get a configuration from ConfigurationList
      ConfigurationPtr getConf(int position) const;

      // returns last configuration 
      ConfigurationPtr back() const;

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
      // TODO should be independent on id sequence
      longint num_of_part;
      longint local_num_of_part;
      longint min_id;
      longint max_id;
      
      string key; // it can be "position", "velocity" or "unfolded"
      
    private:

      void pushConfig(ConfigurationPtr config);
 
      // the list of snapshots
      ConfigurationList configurations;
    };
  }
}

#endif
