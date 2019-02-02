/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)

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
#ifndef _IO_DUMPXTCADRESS_HPP
#define _IO_DUMPXTCADRESS_HPP


#include <boost/serialization/map.hpp>
#include <string>
#include <iostream>
#include "mpi.hpp"

#include "gromacs/fileio/xtcio.h"
#include "types.hpp"
#include "System.hpp"
#include "io/FileBackup.hpp"
#include "ParticleAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "FixedTupleListAdress.hpp"

#include "esutil/Error.hpp"


namespace espressopp {
namespace io {
class DumpXTCAdress : public ParticleAccess {
 public:
  DumpXTCAdress(shared_ptr<System> system,
    shared_ptr<FixedTupleListAdress> _ftpl,
    shared_ptr<integrator::MDIntegrator> _integrator,
    std::string _file_name,
    bool _unfolded,
    real _length_factor,
    bool _append):
      ParticleAccess(system),
      ftpl(_ftpl),
      fio(NULL),
      xtcprec(1000),
      integrator(_integrator),
      file_name(_file_name),
      unfolded(_unfolded),
      length_factor(_length_factor),
      append(_append) {
    if (system->comm->rank() == 0 && !append) {
      FileBackup backup(file_name);  // backup trajectory if it already exists
    }
  }
  ~DumpXTCAdress() { }

  void perform_action() {
    dump();
  }

  void dump();

  std::string getFilename() { return file_name; }
  void setFilename(std::string v) { file_name = v; }
  bool getUnfolded() { return unfolded; }
  void setUnfolded(bool v) { unfolded = v; }
  bool getAppend() { return append; }
  void setAppend(bool v) { append = v; }

  static void registerPython();

 private:
  t_fileio *fio;
  static const int dim = 3;
  real xtcprec;

  // integrator we need to know an integration step
  shared_ptr<integrator::MDIntegrator> integrator;
  shared_ptr<FixedTupleListAdress> ftpl;

  std::string file_name;

  bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
  bool append;  // append to existing trajectory file or create a new one
  real length_factor;

  bool open(const char *mode);
  void close();
  void write(int natoms, int step, float time, Real3D *box, Real3D *x, float prec);
};
}  // end namespace io
}  // end namespace espressopp
#endif
