#include "MeanSquareDispl.hpp"

using namespace std;
//using namespace espresso;

namespace espresso {
  namespace analysis {

    //using namespace iterator;
    
    /*
     * calc <r^2> the output is the average mean sq. displacement over 3 directions.
     * !! Important!! For D calculation one have to divide slope by /2.0 (not 6.0)! It is
     * already taken into account.
     * !! all confs should contain the same number of particles
    */
    
    python::list MeanSquareDispl::compute() const{
      
      int M = getListSize();
      real* totZ;
      totZ = new real[M];
      real* Z;
      Z = new real[M];

      python::list pyli;
      
      System& system = getSystemRef();
      
      // COM calculation
      vector<Real3D> centerOfMassList;
      for(int m=0; m<M; m++){
        Real3D posCOM = Real3D(0.0,0.0,0.0);
        real mass = 0.0;
        Real3D posCOM_sum = Real3D(0.0,0.0,0.0);
        real mass_sum = 0.0;

        for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
          size_t i = itr->first;
          int whichCPU = itr->second;
          
          if(system.comm->rank()==whichCPU){
            Real3D pos = getConf(m)->getCoordinates(i);
            posCOM += pos;
            mass += 1;
          }
        }

        boost::mpi::all_reduce(*mpiWorld, posCOM, posCOM_sum, std::plus<Real3D>());
        boost::mpi::all_reduce(*mpiWorld, mass, mass_sum, std::plus<real>());

        centerOfMassList.push_back( posCOM_sum / mass_sum );
      }
      
      // MSD calculation
      for(int m=0; m<M; m++){
        
        totZ[m] = 0.0;
        Z[m] = 0.0;
        for(int n=0; n<M-m; n++){
          for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
            size_t i = itr->first;
            int whichCPU = itr->second;
            
            if(system.comm->rank()==whichCPU){
              Real3D pos1 = getConf(n + m)->getCoordinates(i) - centerOfMassList[n+m];
              Real3D pos2 = getConf(n)->getCoordinates(i)     - centerOfMassList[n];
              Real3D delta = pos2 - pos1;
              Z[m] += delta.sqr();
            }
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
      
      real coef = 3 * num_of_part;
      
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
