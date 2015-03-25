/*
  Copyright (c) 2015
      Jakub Krajniak (jkrajniak at gmail.com)
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

#include <functional>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "DumpH5MD.hpp"
#include "storage/Storage.hpp"
#include "Version.hpp"
#include "iterator/CellListIterator.hpp"

#include "bc/BC.hpp"
#include "io/ch5md.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

using namespace espressopp::analysis;  //NOLINT

namespace espressopp {
namespace io {

LOG4ESPP_LOGGER(DumpH5MD::theLogger, "DumpH5MD");

DumpH5MD::DumpH5MD(
    shared_ptr<System> system,
    shared_ptr<integrator::MDIntegrator> integrator,
    std::string file_name,
    std::string h5md_group_name,
    bool unfolded,
    std::string author,
    std::string email,
    bool save_force = false,
    bool save_vel = false):
        ParticleAccess(system),
        system_(system),
        integrator_(integrator),
        file_name_(file_name),
        h5md_group_(h5md_group_name),
        unfolded_(unfolded),
        save_force_(save_force),
        save_vel_(save_vel) {
  /// Backups the file.
  if (system->comm->rank() == 0) {
    FileBackup backup(file_name);
  }

  int myN = system_->storage->getNRealParticles();
  int nparticles = 0;
  boost::mpi::all_reduce(*system_->comm, myN, nparticles, std::plus<int>());

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  MPI_Info info = MPI_INFO_NULL;;
  H5Pset_fapl_mpio(plist_id, *system->comm, info);

  Version version;
  file_ = h5md_create_file(file_name.c_str(), author.c_str(), email.c_str(), version.name().c_str(),
      version.version().c_str(), plist_id);
  H5Pclose(plist_id);
  particles_ = h5md_create_particles_group(file_, h5md_group_name.c_str());

  int dims[2];
  dims[0] = nparticles;
  dims[1] = 3;
  particles_.position = h5md_create_time_data(particles_.group, "position", 2, dims,
      H5T_NATIVE_DOUBLE, NULL);
  if (save_vel_) {
    LOG4ESPP_DEBUG(theLogger, "Enabling velocity storage.");
    particles_.velocity = h5md_create_time_data(particles_.group, "velocity", 2, dims,
        H5T_NATIVE_DOUBLE, NULL);
  }
  if (save_force_) {
    LOG4ESPP_DEBUG(theLogger, "Enabling force storage.");
    particles_.force = h5md_create_time_data(particles_.group, "force", 2, dims,
        H5T_NATIVE_DOUBLE, NULL);
  }
  dims[1] = 1;
  particles_.id = h5md_create_time_data(particles_.group, "id", 2, dims, H5T_NATIVE_INT, NULL);
  particles_.species = h5md_create_time_data(particles_.group, "species", 2, dims, H5T_NATIVE_INT,
      NULL);
  dims[1] = 3;

  /// TODO(jakub): Handle also other cases when the box is not periodic.
  const char *boundary[] = {"periodic", "periodic", "periodic"};
  Real3D box = system->bc->getBoxL();
  double edges[3];
  edges[0] = box[0];
  edges[1] = box[1];
  edges[2] = box[2];
  h5md_create_box(&particles_, 3, boundary, false, edges, NULL);

  closed_ = false;
}

void DumpH5MD::dump() {
  /// Get the local particles.
  int myN = system_->storage->getNRealParticles();
  /// Communicate number of particles to calculate offset.
  std::vector<int> all_N(system_->comm->size(), 0);;
  boost::mpi::all_gather(*system_->comm, myN, all_N);

  double *vel = NULL;
  double *force = NULL;
  if (save_force_)
    force = new double[myN*3];
  if (save_vel_)
    vel = new double[myN*3];

  double *r = new double[myN*3];  // buffer for coordinates;
  int *kSp = new int[myN];  // buffer for species;
  int *kIdd = new int[myN];  // buffer for ids;

  /// Fill values with local data;
  CellList realCells = system_->storage->getRealCells();
  int i = 0;
  if (unfolded_) {   // TODO(jakub): Perhaps that could be simplified.
    Real3D box = system_->bc->getBoxL();
    for (iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      kIdd[i] = cit->id();
      kSp[i] = cit->type();
      Real3D &pos = cit->position();
      Int3D &img = cit->image();
      r[i*3] = pos[0] + img[0] * box[0];
      r[i*3 + 1] = pos[1] + img[1] * box[1];
      r[i*3 + 2] = pos[2] + img[2] * box[2];
      if (save_force_) {
        Real3D &f = cit->force();
        force[i*3] = f[0];
        force[i*3 + 1] = f[1];
        force[i*3 + 2] = f[2];
      }
      if (save_vel_) {
        Real3D &v = cit->velocity();
        vel[i*3] = v[0];
        vel[i*3 + 1] = v[1];
        vel[i*3 + 2] = v[2];
      }
      i++;
    }
  } else {
    Real3D box = system_->bc->getBoxL();
    for (iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      kIdd[i] = cit->id();
      kSp[i] = cit->type();
      Real3D &pos = cit->position();
      r[i*3] = pos[0];
      r[i*3 + 1] = pos[1];
      r[i*3 + 2] = pos[2];
      if (save_force_) {
        Real3D &f = cit->force();
        force[i*3] = f[0];
        force[i*3 + 1] = f[1];
        force[i*3 + 2] = f[2];
      }
      if (save_vel_) {
        Real3D &v = cit->velocity();
        vel[i*3] = v[0];
        vel[i*3 + 1] = v[1];
        vel[i*3 + 2] = v[2];
      }
      i++;
    }
  }
  int step = integrator_->getStep();
  double time = integrator_->getTimeStep() * step;

  /// Prepare offset for storing data in parallel
  int offset = 0;
  for (int i_rank = 0; i_rank < system_->comm->rank(); i_rank++)
    offset += all_N[i_rank];

  /// Enabling mpi
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  h5md_append(particles_.position, r, step, time, offset, plist_id, myN);
  h5md_append(particles_.species, kSp, step, time, offset, plist_id, myN);
  h5md_append(particles_.id, kIdd, step, time, offset, plist_id, myN);
  if (save_vel_)
    h5md_append(particles_.velocity, vel, step, time, offset, plist_id, myN);
  if (save_force_)
    h5md_append(particles_.force, force, step, time, offset, plist_id, myN);

  H5Pclose(plist_id);

  delete[] r;
  delete[] kSp;
  delete[] kIdd;
  if (save_force_)
    delete[] force;
  if (save_vel_)
    delete[] vel;
}

void DumpH5MD::close() {
  LOG4ESPP_DEBUG(theLogger, "Closing file.");
  h5md_close_time_data(particles_.position);
  h5md_close_time_data(particles_.id);
  h5md_close_time_data(particles_.species);
  if (save_force_)
    h5md_close_time_data(particles_.force);
  if  (save_vel_)
    h5md_close_time_data(particles_.velocity);
  H5Gclose(particles_.group);
  h5md_close_file(file_);
  closed_ = true;
}

DumpH5MD::~DumpH5MD() {
  LOG4ESPP_DEBUG(theLogger, "Running ~DumpH5MD.");
  if (!closed_)
    close();
}

void DumpH5MD::registerPython() {
  using namespace espressopp::python;  //NOLINT

  class_<DumpH5MD, bases<ParticleAccess>, boost::noncopyable >
  ("io_DumpH5MD", init< shared_ptr< System >,
                       shared_ptr< integrator::MDIntegrator >,
                       std::string,
                       std::string,
                       bool,
                       std::string,
                       std::string,
                       bool,
                       bool>())
    .add_property("file_id", &DumpH5MD::file_id)
    .def("dump", &DumpH5MD::dump)
    .def("close", &DumpH5MD::close);
}
}  // end namespace io
}  // end namespace espressopp
