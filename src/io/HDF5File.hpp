/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & Johannes Gutenberg-Universität Mainz

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
#ifndef _IO_HDF5FILE_HPP
#define _IO_HDF5FILE_HPP

// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "ParticleAccess.hpp"  // keep python.hpp on top
#include "integrator/MDIntegrator.hpp"
#include "io/FileBackup.hpp"
#include "esutil/Error.hpp"
#include <string>

namespace espressopp {
  namespace io{

    class HDF5File : public ParticleAccess {

    public:

      HDF5File(shared_ptr<System> system,
              shared_ptr<integrator::MDIntegrator> _integrator,
              std::string _file_name,
			  int _iomode,
              bool _unfolded,
              real _length_factor,
              std::string _length_unit,
              bool _append) :
                        ParticleAccess(system),
                        integrator(_integrator),
                        file_name( _file_name ),
						iomode(_iomode),
                        unfolded(_unfolded),
                        length_factor(_length_factor),
                        append(_append){
        setLengthUnit(_length_unit);

        if (iomode == 1 || iomode == 0) {

			if (system->comm->rank() == 0  && !append) {
				FileBackup backup(file_name);
			}
        } else if (iomode == 2) {

        	if (!append) {
				int rank = system->comm->rank();
				size_t filename_length = file_name.length();
				std::string suffix = file_name.substr(filename_length-3, 3);
				std::string base_filename = file_name.substr(0,filename_length-3);
				std::string rankstring = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) )->str();
				std::string final_name = base_filename + "_" + rankstring + suffix;
				FileBackup backup(final_name);
        	}

        }
      }
      ~HDF5File() {}

      void perform_action(){
        write();
      }

      void write_n_to_1();
      void write_n_to_n();

      //void read_n_to_1();
      //void read_n_to_n();

      void write();
//      void read() {
//    	  if (iomode == 1 || iomode == 0) {
//    		  read_n_to_1();
//    	  } else if (iomode == 2) {
//    		  read_n_to_n();
//    	  }
//      }

      std::string getFilename(){return file_name;}
      void setFilename(std::string v){file_name = v;}
      inline int getIomode(){return iomode;}
      inline void setIomode(int v){iomode = v;}
      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}
      bool getAppend(){return append;}
      void setAppend(bool v){append = v;}

      std::string getLengthUnit(){return length_unit;}
      void setLengthUnit(std::string v){
        esutil::Error err( getSystem()->comm );
        if( v != "LJ" && v != "nm" && v != "A" ){
          std::stringstream msg;
          msg<<"Wrong unit length: "<< v << "  It should be string: LJ, nm or A" <<"\n";
          err.setException( msg.str() );
          err.checkException();
        }

        length_unit = v;
      }
      real getLengthFactor(){return length_factor;}
      void setLengthFactor(real v){length_factor = v;}

      static void registerPython();

    protected:

      //static LOG4ESPP_DECL_LOGGER(logger);

    private:

      // integrator we need to know an integration step
      shared_ptr<integrator::MDIntegrator> integrator;

      std::string file_name;
      int iomode; // 0: serial, 1: N-to-1, 2: N-to-N; real 0 not there now

      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
      bool append; //append to existing trajectory file or create a new one
      real length_factor;  // for example
      std::string length_unit; // length unit: {could be LJ, nm, A} it is just for user info
    };
  }
}

#endif

