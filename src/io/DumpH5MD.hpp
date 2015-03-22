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

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

#include "esutil/Error.hpp"

#include "io/ch5md.hpp"

using namespace espressopp::analysis;  // NOLINT

namespace espressopp {
namespace io {

class DumpH5MD : public ParticleAccess {
 public:
  DumpH5MD(shared_ptr<System> system,
           shared_ptr<integrator::MDIntegrator> integrator,
           std::string file_name,
           std::string h5md_group_name,
           bool unfolded,
           std::string author);

  ~DumpH5MD();

  void perform_action() { dump(); }
  hid_t file_id() { return file_.id; }
  void dump();
  void close();
  static void registerPython();

 private:
  // integrator we need to know an integration step
  shared_ptr<integrator::MDIntegrator> integrator_;
  shared_ptr<System> system_;

  std::string file_name_;
  std::string h5md_group_;

  /// If set to true the the coordinates are unfolded. By default it is folded.
  bool unfolded_;

  /// Name of the author of file.
  std::string author_;

  /// Number of particles;
  int nparticles_;

  /// H5MD
  h5md_file file_;
  h5md_particles_group particles_;
};

}  // end namespace io
}  // end namespace espressopp
#endif
