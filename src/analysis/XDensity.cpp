#include "python.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "Configuration.hpp"
#include "XDensity.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"
#include "math.h"

#include <boost/serialization/map.hpp>

using namespace espresso;
using namespace iterator;
using namespace std;

namespace espresso {
  namespace analysis {
    
    // splitN is a level of discretisation of density profile (how many elements it contains)
    python::list XDensity::computeArray(int splitN) const {

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      Real3D Li = system.bc->getBoxL();
      
      int nprocs = system.comm->size();
      int myrank = system.comm->rank();
      
      real histogram [splitN];
      for(int i=0;i<splitN;i++) histogram[i]=0;
          
      real dr = Li[0] / (real)splitN;
          
      int num_part = 0;
      ConfigurationPtr config = make_shared<Configuration> ();
      for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, Real3D > conf;
        if (rank_i == myrank) {
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            int id = cit->id();
            conf[id] = cit->position();
          }
    	}

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        // for simplicity we will number the particles from 0
        for (map<size_t,Real3D>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
          //size_t id = itr->first;
          Real3D p = itr->second;
          config->set(num_part, p[0], p[1], p[2]);
          num_part ++;
        }
      }
      // now all CPUs have all particle coords and num_part is the total number of particles
      
      // use all cpus
      // TODO it could be a problem if   n_nodes > num_part
      int numi = num_part / nprocs + 1;
      int mini = myrank * numi;
      int maxi = mini + numi;
      
      if(maxi>num_part) maxi = num_part;
      
      int perc=0, perc1=0;
      for(int i = mini; i<maxi; i++){
        Real3D coordP1 = config->getCoordinates(i);
        int bin = floor (coordP1[0]/dr);
        histogram[bin] += 1.0;
        
        /*if(system.comm->rank()==0){
          perc = (int)(100*(real)(i-mini)/(real)(maxi-mini));
          if(perc>perc1){
            cout<<"calculation progress (density profile along x axis): "<< perc << " %\r"<<flush;
            perc1 = perc;
          }
        }*/
      }
      //if(system.comm->rank()==0)
        //cout<<"calculation progress (density profile along x axis): "<< 100 << " %" <<endl;

      real totHistogram[splitN];
      boost::mpi::all_reduce(*system.comm, histogram, splitN, totHistogram, plus<real>());
      
      // normalizing
      int nconfigs = 1; //config - 1
      real rho = (real)num_part / splitN;
      
      for(int i=0; i < splitN; i++){
        totHistogram[i] /= rho;
      }
      
      python::list pyli;
      for(int i=0; i < splitN; i++){
        pyli.append( totHistogram[i] );
      }
      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real XDensity::compute() const {
      return -1.0;
    }

    void XDensity::registerPython() {
      using namespace espresso::python;
      class_<XDensity, bases< Observable > >
        ("analysis_XDensity", init< shared_ptr< System > >())
        .def("compute", &XDensity::computeArray)
      ;
    }
  }
}
