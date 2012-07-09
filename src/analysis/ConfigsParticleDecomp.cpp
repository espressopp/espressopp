#include "mpi.h"

#include "ConfigsParticleDecomp.hpp"
#include "iterator/CellListIterator.hpp"

#include <boost/serialization/map.hpp>

using namespace std;
using namespace espresso;

namespace espresso {
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
        cout << "Error. Velocities::get <out-of-range>" << endl;
        return shared_ptr<Configuration>();
      }
    }

    ConfigurationPtr ConfigsParticleDecomp::back() const{
      return configurations.back();
    }

    void ConfigsParticleDecomp::pushConfig(ConfigurationPtr config){
      configurations.push_back(config);
    }

    void ConfigsParticleDecomp::gather() {
      System& system = getSystemRef();
      
      int nprocs = system.comm->size();
      int myrank = system.comm->rank();

      if(myrank==0){
        // check whether the number of particles is the same during the gathering
        int curNumP = system.storage->getNRealParticles();
        if(curNumP != num_of_part){
          cout<<"   ConfigsParticleDecomp gathers the configurations of the same system\n"
                " with the same number of particles. If you need to store the systems\n"
                " with different number of particles you should use something else."
                " E.g `Configurations`"<< endl;
          return;
        }
      }

      ConfigurationPtr config = make_shared<Configuration> ();
      for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, Real3D > conf;
        if (rank_i == myrank) {
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            int id = cit->id();
            Real3D& vel = cit->velocity();
            conf[id] = vel;
          }
    	}

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        for (map<size_t,Real3D>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
          size_t id = itr->first;
          Real3D p = itr->second;
          if(id>=min_id && id<max_id)
            config->set(id, p[0], p[1], p[2]);
        }
      }

      pushConfig(config);
    }
    
    // Python wrapping
    void ConfigsParticleDecomp::registerPython() {
      using namespace espresso::python;

      class_<ConfigsParticleDecomp, boost::noncopyable >(
        "analysis_ConfigsParticleDecomp", no_init
        //init< shared_ptr< System > >()
      )
      .def_readonly("size", &ConfigsParticleDecomp::getListSize)
      
      .def("gather", &ConfigsParticleDecomp::gather)
      .def("__getitem__", &ConfigsParticleDecomp::getConf)
      .def("back", &ConfigsParticleDecomp::back)
      .def("all", &ConfigsParticleDecomp::all)
      .def("clear", &ConfigsParticleDecomp::clear)
      .def("compute", &ConfigsParticleDecomp::compute)
      ;
    }
  }
}
