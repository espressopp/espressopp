#include "mpi.h"

#include "Autocorrelation.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"

using namespace std;

namespace espresso {
  namespace analysis {

    using namespace iterator;

    int Autocorrelation::getListSize() const{
      return valueList.size();
    }
   
    vector<Real3D> Autocorrelation::all() const{
      return valueList;
    }

    Real3D Autocorrelation::getValue(int position) const{
      int nconfigs = getListSize();
      if (0 <= position and position < nconfigs) {
        return valueList[position];
      }
      else{
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        stringstream msg;
        msg << "Error. Velocities::get <out-of-range>";
        err.setException( msg.str() );
        return Real3D(0,0,0);
      }
    }

    void Autocorrelation::pushValue(Real3D value){
      valueList.push_back(value);
    }

    void Autocorrelation::gather(Real3D value) {
      pushValue(value);
    }
    
    python::list Autocorrelation::compute() {
      int M = getListSize();
      
      python::list pyli;
      
      System& system = getSystemRef();
      
      int n_nodes = system.comm -> size();
      int this_node = system.comm -> rank();

      // TODO it could be a problem if   n_nodes > total_num !!!
      int num_m[n_nodes];
      int num_mH[n_nodes];
      
      if(this_node == 0){
        
        // it is 1+2+3+...+M
        double local_num = ( (double)M * (double)(M + 1) / 2.0) / (double)n_nodes + 1.0;
        
        for(int i=0; i<n_nodes; i++){
          double max_num = (i+1) * local_num;
          
          int m_max = (int) ( ( sqrt(1.0+8.0*max_num) - 1.0 )/2. );
          
          int lastNum = (i==0) ? 0 : num_mH[i-1];
          
          if(m_max-lastNum==0)
            num_mH[i] = num_mH[i-1] + 1;
          else
            num_mH[i] = m_max;
          
          if(num_mH[i] > M) num_mH[i] = M;
        }
        
        for(int i=0; i<n_nodes-1; i++)
          num_m[i] = M - num_mH[n_nodes-2-i];
        num_m[n_nodes-1] = M;
      }

      boost::mpi::broadcast(*system.comm, num_m, n_nodes, 0);
      
      // now num_m[i], i - cpu number, is a number of series in "for" statement for each cpu
      
      int min_m = (this_node==0) ? 0 : num_m[this_node-1];
      int max_m = num_m[this_node];

      real* Z;
      Z = new real[M];
      
      for(int m=min_m; m<max_m; m++){
        Z[m] = 0.0;
        for(int n=0; n<M-m; n++){
          Z[m] += getValue(n + m) *  getValue(n);
        }
        if(system.comm->rank()==0)
          cout<<"calculation progress (autocorrelation): "<< (int)(100*(real)m/(real)(max_m-min_m)) << " %\r" <<flush;
      }
      
      if(system.comm->rank()==0)
        cout<<"calculation progress (autocorrelation): 100 %" <<endl;
      
      real coef = 3.0; // only if value is Real3D
      
      for(int m=min_m; m<max_m; m++){
        Z[m] /= ( (real)(M-m)*coef );
      }
      
      // TODO probably could be done nicer. gather doesn't work with different length of array
      unsigned long int MM = M * n_nodes;
      real* totZ = new real[MM];
      boost::mpi::gather(*system.comm, Z, M, totZ, 0);
      
      if(this_node == 0){
        int count = 0;
        for(int m=0; m<M; m++){
          if( m >= num_m[count] ) count ++;
          pyli.append( totZ[M*count + m] );
        }
      }
      
      delete [] Z;
      Z = NULL;
      delete [] totZ;
      totZ = NULL;
      
      return pyli;
    }
    
    // Python wrapping
    void Autocorrelation::registerPython() {
      using namespace espresso::python;

      class_<Autocorrelation>(
        "analysis_Autocorrelation",
        init< shared_ptr< System > >()
      )
      .def_readonly("size", &Autocorrelation::getListSize)
      
      .def("gather", &Autocorrelation::gather)
      .def("__getitem__", &Autocorrelation::getValue)
      .def("all", &Autocorrelation::all)
      .def("clear", &Autocorrelation::clear)
      .def("compute", &Autocorrelation::compute)
      
      ;
    }
  }
}
