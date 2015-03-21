/*
  Copyright (C) 2012,2013
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

// ESPP_CLASS
#ifndef _IO_DumpH5MD_HPP
#define _IO_DumpH5MD_HPP

#include <hdf5.h>
#include <boost/serialization/map.hpp>
#include <string>
#include "mpi.hpp"
#include "types.hpp"
#include "System.hpp"
#include "io/FileBackup.hpp"
#include "ParticleAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

#include "esutil/Error.hpp"


namespace espresso {
namespace io {

class DumpH5MD : public ParticleAccess {
 public:
  DumpH5MD(shared_ptr<System> system,
           shared_ptr<integrator::MDIntegrator> integrator,
           std::string file_name,
           std::string h5md_group_name,
           bool unfolded,
           std::string length_unit
           std::string author):
               ParticleAccess(system),
               integrator_(integrator),
               file_name_(file_name),
               h5md_group_(h5md_group_name),
               unfolded_(unfolded),
               length_unit_(length_unit),
               author_(author) {
    if (system->comm->rank() == 0) {
      FileBackup backup(file_name);
    }

    /** Prepare the file. */
    file_id_ = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  ~DumpH5MD() { }

  void perform_action() { dump(); }
  void dump();
  static void registerPython();

 private:
  // integrator we need to know an integration step
  shared_ptr<integrator::MDIntegrator> integrator;

  std::string file_name_;
  std::string h5md_group_;

  /// If set to true the the coordinates are unfolded. By default it is folded.
  bool unfolded_;
  /// Attributes with length units.
  std::string length_unit_;

  /// Name of the author of file.
  std::string author_;

  /// H5MD
  hid_t file_id_;
  hid_t particles_;
};

}  // end namespace io
}  // end namespace espresso
#endif
