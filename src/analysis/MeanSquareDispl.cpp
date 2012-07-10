#include "MeanSquareDispl.hpp"

using namespace std;
//using namespace espresso;

namespace espresso {
  namespace analysis {

    //using namespace iterator;
    
    // calc <vx(0) * vx(t)>
    // !!!! all confs should contain the same num  of particles
    python::list MeanSquareDispl::compute() const{
      
      int M = getListSize();
      real* totZ;
      totZ = new real[M];
      real* Z;
      Z = new real[M];

      python::list pyli;
      
      System& system = getSystemRef();
      
      for(int m=0; m<M; m++){
        totZ[m] = 0.0;
        Z[m] = 0.0;
        for(int n=0; n<M-m; n++){
          for(int i=min_id; i<max_id; i++){
            Real3D pos1 = getConf(n + m)->getCoordinates(i);
            Real3D pos2 = getConf(n)->getCoordinates(i);
            Real3D delta = pos2 - pos1;
            Z[m] += delta.sqr();
          }
        }
        if(system.comm->rank()==0)
          cout<<"calculation progress (mean square displacement): "<< (int)(100*(real)m/(real)M) << " %\r"<<flush;
      }
      if(system.comm->rank()==0)
        cout<<"calculation progress (mean square displacement): "<< (int)(100*(real)M/(real)M) << " %" <<endl;
      
      boost::mpi::all_reduce( *system.comm, Z, M, totZ, plus<real>() );
      
      for(int m=0; m<M; m++){
        totZ[m] /= (real)(M - m);
      }
      
      real coef =  num_of_part * 3;
      
      for(int m=0; m<M; m++){
        totZ[m] /= coef;
        pyli.append( totZ[m] );
      }
      
      delete [] Z;
      Z = NULL;
      delete [] totZ;
      totZ = NULL;
      
      return pyli;
    }
    
    // Python wrapping

    void MeanSquareDispl::registerPython() {
      using namespace espresso::python;

      class_<MeanSquareDispl, bases<ConfigsParticleDecomp> >(
        "analysis_MeanSquareDispl",
        init< shared_ptr< System > >()
      );
    }
  }
}
