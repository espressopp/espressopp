/*
  Copyright (C) 2021
      Sebastian Eibl, Max Planck Computing & Data Facility

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

#include "RestoreH5MDParallel.hpp"

#include "iterator/CellListIterator.hpp"
#include "storage/Storage.hpp"

namespace espressopp
{
namespace io
{
template <typename T>
void RestoreH5MDParallel::readParallel(hid_t fileId,
                                       const std::string& dataset,
                                       std::vector<T>& data)
{
    auto dset = CHECK_HDF5(H5Dopen(fileId, dataset.c_str(), H5P_DEFAULT));
    auto dspace = CHECK_HDF5(H5Dget_space(dset));

    // get global dimensions
    std::vector<hsize_t> globalDims;
    auto ndims = CHECK_HDF5(H5Sget_simple_extent_ndims(dspace));
    CHECK_GREATER(ndims, 0);
    globalDims.resize(ndims);
    CHECK_HDF5(H5Sget_simple_extent_dims(dspace, globalDims.data(), nullptr));

    // get local dimensions and offset
    std::vector<hsize_t> localDims = globalDims;
    hsize_t localOffset = 0;
    for (auto rk = 0; rk < rank; ++rk)
    {
        localOffset += globalDims[1] / uint_c(numProcesses) +
                       (globalDims[1] % uint_c(numProcesses) > uint_c(rk) ? 1ul : 0ul);
    }
    auto localSize = globalDims[1] / uint_c(numProcesses) +
                     (globalDims[1] % uint_c(numProcesses) > uint_c(rank) ? 1ul : 0ul);
    localDims[0] = 1;  // only read one timeframe
    localDims[1] = localSize;

    // set up local part of the input file
    std::vector<hsize_t> offset(globalDims.size(), 0);
    offset[1] = localOffset;
    std::vector<hsize_t> stride(globalDims.size(), 1);
    std::vector<hsize_t> count(globalDims.size(), 1);
    // check if in bounds
    for (auto i = 0; i < int_c(globalDims.size()); ++i)
    {
        CHECK_LESS_EQUAL(localDims[i] + offset[i], globalDims[i], "i = " << i);
    }
    auto fileSpace = CHECK_HDF5(H5Dget_space(dset));
    CHECK_HDF5(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), stride.data(),
                                   count.data(), localDims.data()));

    // set up memory data layout
    hid_t memSpace =
        CHECK_HDF5(H5Screate_simple(int_c(localDims.size()), localDims.data(), nullptr));
    auto linLocalSize =
        std::accumulate(localDims.begin(), localDims.end(), hsize_t(1), std::multiplies<>());
    data.resize(linLocalSize);

    // read
    auto dataread = CHECK_HDF5(H5Pcreate(H5P_DATASET_XFER));
    CHECK_HDF5(H5Pset_dxpl_mpio(dataread, H5FD_MPIO_COLLECTIVE));
    CHECK_HDF5(H5Dread(dset, typeToHDF5<T>(), memSpace, fileSpace, dataread, data.data()));

    // close
    CHECK_HDF5(H5Sclose(fileSpace));
    CHECK_HDF5(H5Sclose(memSpace));
    CHECK_HDF5(H5Pclose(dataread));
    CHECK_HDF5(H5Sclose(dspace));
    CHECK_HDF5(H5Dclose(dset));
}

void RestoreH5MDParallel::updateCache()
{
    boost::mpi::communicator world;
    comm = world;
    numProcesses = world.size();
    rank = world.rank();
}

void RestoreH5MDParallel::restore()
{
    updateCache();

    auto info = MPI_INFO_NULL;

    auto plist = CHECK_HDF5(H5Pcreate(H5P_FILE_ACCESS));
    CHECK_HDF5(H5Pset_fapl_mpio(plist, comm, info));

    auto fileId = CHECK_HDF5(H5Fopen(filename_.c_str(), H5F_ACC_RDONLY, plist));
    std::vector<int64_t> id;
    if (restoreId)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + idDataset + "/value", id);
        CHECK_EQUAL(id.size() * 1, id.size());
    }
    std::vector<int64_t> type;
    if (restoreType)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + typeDataset + "/value",
                     type);
        CHECK_EQUAL(id.size() * 1, type.size());
    }
    std::vector<double> mass;
    if (restoreMass)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + massDataset + "/value",
                     mass);
        CHECK_EQUAL(id.size() * 1, mass.size());
    }
    std::vector<double> q;
    if (restoreQ)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + qDataset + "/value", q);
        CHECK_EQUAL(id.size() * 1, q.size());
    }
    std::vector<int8_t> ghost;
    if (restoreGhost)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + ghostDataset + "/value",
                     ghost);
        CHECK_EQUAL(id.size() * 1, ghost.size());
    }
    std::vector<double> position;
    if (restorePosition)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + positionDataset + "/value",
                     position);
        CHECK_EQUAL(id.size() * 3, position.size());
    }
    std::vector<double> velocity;
    if (restoreVelocity)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + velocityDataset + "/value",
                     velocity);
        CHECK_EQUAL(id.size() * 3, velocity.size());
    }
    std::vector<double> force;
    if (restoreForce)
    {
        readParallel(fileId, "/particles/" + particleGroupName + "/" + forceDataset + "/value",
                     force);
        CHECK_EQUAL(id.size() * 3, force.size());
    }

    for (auto i = 0; i < int_c(id.size()); ++i)
    {
        auto pos = Real3D(position[i * 3 + 0], position[i * 3 + 1], position[i * 3 + 2]);
        auto particle = system_->storage->addParticle(int_c(id[i]), pos, false);
        CHECK_NOT_NULLPTR(particle, "particle creation was rejected!");
        if (restoreId) particle->id() = id[i];
        if (restoreType) particle->type() = type[i];
        if (restoreMass) particle->mass() = mass[i];
        if (restoreQ) particle->q() = q[i];
        if (restoreGhost) particle->ghost() = ghost[i];
        if (restorePosition)
        {
            particle->position()[0] = position[i * 3 + 0];
            particle->position()[1] = position[i * 3 + 1];
            particle->position()[2] = position[i * 3 + 2];
        }
        if (restoreVelocity)
        {
            particle->velocity()[0] = velocity[i * 3 + 0];
            particle->velocity()[1] = velocity[i * 3 + 1];
            particle->velocity()[2] = velocity[i * 3 + 2];
        }
        if (restoreForce)
        {
            particle->force()[0] = force[i * 3 + 0];
            particle->force()[1] = force[i * 3 + 1];
            particle->force()[2] = force[i * 3 + 2];
        }
    }

    CHECK_HDF5(H5Fclose(fileId));
}

void RestoreH5MDParallel::registerPython()
{
    using namespace espressopp::python;

    class_<RestoreH5MDParallel>("io_RestoreH5MDParallel", init<shared_ptr<System>, std::string>())
        .def_readwrite("particleGroupName", &RestoreH5MDParallel::particleGroupName)
        .def_readwrite("restoreId", &RestoreH5MDParallel::restoreId)
        .def_readwrite("restoreType", &RestoreH5MDParallel::restoreType)
        .def_readwrite("restoreMass", &RestoreH5MDParallel::restoreMass)
        .def_readwrite("restoreQ", &RestoreH5MDParallel::restoreQ)
        .def_readwrite("restoreGhost", &RestoreH5MDParallel::restoreGhost)
        .def_readwrite("restorePosition", &RestoreH5MDParallel::restorePosition)
        .def_readwrite("restoreVelocity", &RestoreH5MDParallel::restoreVelocity)
        .def_readwrite("restoreForce", &RestoreH5MDParallel::restoreForce)
        .def_readwrite("idDataset", &RestoreH5MDParallel::idDataset)
        .def_readwrite("typeDataset", &RestoreH5MDParallel::typeDataset)
        .def_readwrite("massDataset", &RestoreH5MDParallel::massDataset)
        .def_readwrite("qDataset", &RestoreH5MDParallel::qDataset)
        .def_readwrite("ghostDataset", &RestoreH5MDParallel::ghostDataset)
        .def_readwrite("positionDataset", &RestoreH5MDParallel::positionDataset)
        .def_readwrite("velocityDataset", &RestoreH5MDParallel::velocityDataset)
        .def_readwrite("forceDataset", &RestoreH5MDParallel::forceDataset)
        .def("restore", &RestoreH5MDParallel::restore);
}
}  // namespace io
}  // namespace espressopp