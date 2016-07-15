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
#include <fstream>
#include <iomanip>
#include "DumpGROAdress.hpp"
#include "storage/Storage.hpp"

#include "bc/BC.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"
#include "analysis/ConfigurationsExtAdress.hpp"

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;

namespace espressopp {
  namespace io {
      
    void DumpGROAdress::dump(){
      
      shared_ptr<System> system = getSystem();
      ConfigurationsExtAdress conf( system,fixedTupleList );
      conf.setUnfolded(unfolded);
      conf.gather();
      
      if( system->comm->rank()==0 ){
        ConfigurationExtPtr conf_real = conf.back();
        
        int num_of_particles = conf_real->getSize();
        //int dimension=0;
        //if (num_of_particles != 0 )
          //dimension = 6;//conf_real->getProperties(0).getDimension();
      
        char *ch_f_name = new char[file_name.length() + 1];
        strcpy(ch_f_name, file_name.c_str());
        ofstream myfile (ch_f_name, ios::out | ios::app);
        if (myfile.is_open()){
          //GRO file format, see http://manual.gromacs.org/online/gro.html
          //first line: system description
          //second line: number of particles
          //line 3-n+2: particles
          //line n+3: box info
          //repeat for each frame

          //myfile << num_of_particles << endl;
          myfile << setiosflags(ios::fixed); // needed for fixed-width output
          myfile << "system description, "
                 << "current step="<< integrator->getStep() << ", "
                 << "length unit=" << length_unit << endl;
          myfile << setw(5) << num_of_particles << endl;
          
          
          // for noncubic simulation boxes
          //myfile << Li[0] * length_factor << "  0.0  0.0  0.0  "<< 
           //       Li[1] * length_factor << "  0.0  0.0  0.0  "<< Li[2] * length_factor;
          // additional info to comment line
          //myfile << "  currentStep " << integrator->getStep() << "  lengthUnit "<< length_unit << endl;
         
          //do I need the if statement for length_factor?

          ConfigurationExtIterator cei = conf_real-> getIterator();
          RealND::iterator ii;
          short ind;
          size_t atomnumber;
          for(size_t i=0; i<num_of_particles; i++){
          // "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
            //myfile << setw(5) << i+1;    //FIXME this should be the molecule number, not atom number
            myfile << setw(5) << 10000;
            //myfile << setiosflags(ios::left) << setw(1) << "T" << setw(4) <<   particleIDToType.find(i+1)->second <<resetiosflags(ios::left);  // pid starts at 1 // set(1)+set(4) makes in total 5, as required by the fixed format, should be resname not atomtype
            //stringstream ss;
            //ss << particleIDToType.find(i+1)->second;
            //myfile << setiosflags(ios::right) << setw(5) << (string("T") + ss.str()) << resetiosflags(ios::right);
            myfile << setiosflags(ios::left) << setw(5) << "TTT" <<resetiosflags(ios::left);  // pid starts at 1 
            myfile << setiosflags(ios::right) << setw(5) << (string("TTT")) << resetiosflags(ios::right);
            if ((i+1)<100000) {
              atomnumber=i+1;
            } else {
              atomnumber=i-99999;
            }
            myfile << setw(5) << atomnumber;    //NOTE this is the actual atom number - wrapped at 99999
            // while get token 
            // print with setw(8) << setprecision(3)
            // if more than 3
            // change precision to 4
            RealND line(cei.nextProperties()); //FIXME create line every atom?
            ii = line.begin();
            ind=0; //FIXME ugly! how do I know the number of a loop when using iterators?
            while (ii != line.end() && ind < 3){
              myfile << setw(8) << setprecision(3) << length_factor * *ii;
              ii++;
              ind++;
            }
            while (ii != line.end()){
              myfile << setw(8) << setprecision(4) << length_factor * *ii;
              ii++;
            }
            myfile << endl;

          }
          Real3D Li = system->bc->getBoxL();
          myfile << setw(10) << setprecision(5) << Li[0] * length_factor 
                 << setw(10) << setprecision(5) << Li[1] * length_factor
                 << setw(10) << setprecision(5)  << Li[2] * length_factor
                 << endl;
          myfile.close();
        }
        else cout << "Unable to open file: "<< file_name <<endl;

        delete [] ch_f_name;
      }
    }
      
    // Python wrapping
    void DumpGROAdress::registerPython() {

      using namespace espressopp::python;

      class_<DumpGROAdress, bases<ParticleAccess>, boost::noncopyable >
      ("io_DumpGROAdress", init< shared_ptr< System >, shared_ptr<FixedTupleListAdress>,
                           shared_ptr< integrator::MDIntegrator >, 
                           std::string, 
                           bool,
                           real,
                           std::string,
                           bool>())
        .add_property("filename", &DumpGROAdress::getFilename, 
                                  &DumpGROAdress::setFilename)
        .add_property("unfolded", &DumpGROAdress::getUnfolded, 
                                  &DumpGROAdress::setUnfolded)
        .add_property("append", &DumpGROAdress::getAppend, 
                                  &DumpGROAdress::setAppend)
        .add_property("length_factor", &DumpGROAdress::getLengthFactor, 
                                       &DumpGROAdress::setLengthFactor)
        .add_property("length_unit", &DumpGROAdress::getLengthUnit, 
                                     &DumpGROAdress::setLengthUnit)
        .def("dump", &DumpGROAdress::dump)
      ;
    }
  }
}
