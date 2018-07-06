/*
  Copyright (C) 2017
      Gregor Deichmann (TU Darmstadt, deichmann(at)cpc.tu-darmstadt.de) 
 
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
#ifndef _IO_DUMPXTC_HPP
#define _IO_DUMPXTC_HPP

#include "mpi.hpp"
#include <boost/serialization/map.hpp>
#include "types.hpp"
#include "System.hpp"
#include "io/FileBackup.hpp"
#include "ParticleAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

#include "esutil/Error.hpp"

#include <string>
#include <iostream>

#include "gromacs/fileio/xtcio.h"

namespace espressopp {
  namespace io{

    class DumpXTC : public ParticleAccess {

    public:

      DumpXTC(shared_ptr<System> system, 
              shared_ptr<integrator::MDIntegrator> _integrator,
              std::string _file_name,
              bool _unfolded,
              real _length_factor,
              bool _append):
                        ParticleAccess(system),
                        fio(NULL),
                        xtcprec(1000),
                        integrator(_integrator),
                        file_name( _file_name ),
                        unfolded(_unfolded),
                        length_factor(_length_factor),
                        append(_append){

        if( system->comm->rank()==0 && !append){
          FileBackup backup(file_name); //backup trajectory if it already exists
        }

      }
      ~DumpXTC() {std::cout << "DumpXTC destructor" << std::endl;} // never called, right?

      void perform_action(){
        dump();
      }
      
      void dump();
      
      std::string getFilename(){return file_name;}
      void setFilename(std::string v){file_name = v;}
      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}
      bool getAppend(){return append;}
      void setAppend(bool v){append = v;}

      static void registerPython();
    
    private:
      
      t_fileio *fio;
      static const int dim=3;
      real xtcprec;

      // integrator we need to know an integration step
      shared_ptr<integrator::MDIntegrator> integrator;
      
      std::string file_name;

      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
      bool append; //append to existing trajectory file or create a new one
      real length_factor;   

      bool open(const char *mode);
      void close();
      void write(int natoms,int step,float time,Real3D *box,Real3D *x,float prec);

    };
  }
}

#endif
