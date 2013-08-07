#include <fstream>

#include "DumpXYZ.hpp"
#include "storage/Storage.hpp"

#include "bc/BC.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

using namespace espresso;
using namespace espresso::analysis;
using namespace std;

namespace espresso {
  namespace io {
      
    void DumpXYZ::dump(){
      
      shared_ptr<System> system = getSystem();
      ConfigurationsExt conf( system );
      conf.setUnfolded(unfolded);
      conf.gather();
      
      if( system->comm->rank()==0 ){
        ConfigurationExtPtr conf_real = conf.back();
        
        int num_of_particles = conf_real->getSize();
      
        char *ch_f_name = new char[file_name.length() + 1];
        strcpy(ch_f_name, file_name.c_str());
        ofstream myfile (ch_f_name, ios::out | ios::app);
        if (myfile.is_open()){
          myfile << num_of_particles << endl;
          
          Real3D Li = system->bc->getBoxL();
          
          // for noncubic simulation boxes
          myfile << Li[0] << "  0.0  0.0  0.0  "<< Li[1] << "  0.0  0.0  0.0  "<< Li[2] << endl;
          
          ConfigurationExtIterator cei = conf_real-> getIterator();
          for(size_t i=0; i<num_of_particles; i++){
            myfile << "  0  " << cei.nextCoordinates() << endl;
          }
          myfile.close();
        }
        else cout << "Unable to open file: "<< file_name <<endl;

        delete [] ch_f_name;
      }
    }
      
    // Python wrapping
    void DumpXYZ::registerPython() {

      using namespace espresso::python;

      class_<DumpXYZ, bases<ParticleAccess>, boost::noncopyable >
      ("io_DumpXYZ", init< shared_ptr< System >, 
                           shared_ptr< integrator::MDIntegrator >, 
                           std::string, bool >())
        .add_property("filename", &DumpXYZ::getFilename, 
                                  &DumpXYZ::setFilename)
        .add_property("unfolded", &DumpXYZ::getUnfolded, 
                                  &DumpXYZ::setUnfolded)
        .def("dump", &DumpXYZ::dump)
      ;
    }
  }
}
