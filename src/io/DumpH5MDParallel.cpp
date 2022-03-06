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

#include "DumpH5MDParallel.hpp"

#include "bc/BC.hpp"
#include "iterator/CellListIterator.hpp"
#include "storage/Storage.hpp"
#include "Version.hpp"

namespace espressopp
{
namespace io
{
template <typename T>
void DumpH5MDParallel::writeParallel(hid_t fileId,
                                     const std::string& name,
                                     const std::vector<hsize_t>& globalDims,
                                     const std::vector<hsize_t>& localDims,
                                     const std::vector<T>& data)
{
    CHECK_EQUAL(globalDims.size(), localDims.size());
    CHECK_EQUAL(data.size(), std::accumulate(localDims.begin(), localDims.end(), hsize_t(1),
                                             std::multiplies<>()));

    auto dataspace =
        CHECK_HDF5(H5Screate_simple(int_c(globalDims.size()), globalDims.data(), nullptr));
    auto dataset = CHECK_HDF5(H5Dcreate(fileId, name.c_str(), typeToHDF5<T>(), dataspace,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    std::vector<hsize_t> offset(globalDims.size(), 0);
    offset[1] = particleOffset;
    std::vector<hsize_t> stride(globalDims.size(), 1);
    std::vector<hsize_t> count(globalDims.size(), 1);
    for (auto i = 0; i < int_c(globalDims.size()); ++i)
    {
        CHECK_LESS_EQUAL(localDims[i] + offset[i], globalDims[i], "i = " << i);
    }
    auto dstSpace = CHECK_HDF5(H5Dget_space(dataset));
    CHECK_HDF5(H5Sselect_hyperslab(dstSpace, H5S_SELECT_SET, offset.data(), stride.data(),
                                   count.data(), localDims.data()));

    std::vector<hsize_t> localOffset(globalDims.size(), 0);
    auto srcSpace =
        CHECK_HDF5(H5Screate_simple(int_c(localDims.size()), localDims.data(), nullptr));
    CHECK_HDF5(H5Sselect_hyperslab(srcSpace, H5S_SELECT_SET, localOffset.data(), stride.data(),
                                   count.data(), localDims.data()));

    auto datawrite = CHECK_HDF5(H5Pcreate(H5P_DATASET_XFER));
    CHECK_HDF5(H5Pset_dxpl_mpio(datawrite, H5FD_MPIO_COLLECTIVE));
    CHECK_HDF5(H5Dwrite(dataset, typeToHDF5<T>(), srcSpace, dstSpace, datawrite, data.data()));

    CHECK_HDF5(H5Pclose(datawrite));
    CHECK_HDF5(H5Sclose(dstSpace));
    CHECK_HDF5(H5Sclose(srcSpace));
    CHECK_HDF5(H5Dclose(dataset));
    CHECK_HDF5(H5Sclose(dataspace));
}

void DumpH5MDParallel::writeHeader(hid_t fileId) const
{
    auto group1 = CHECK_HDF5(H5Gcreate(fileId, "/h5md", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    auto group2 =
        CHECK_HDF5(H5Gcreate(fileId, "/h5md/author", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    auto group3 =
        CHECK_HDF5(H5Gcreate(fileId, "/h5md/creator", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    CHECK_HDF5(H5Gclose(group1));
    CHECK_HDF5(H5Gclose(group2));
    CHECK_HDF5(H5Gclose(group3));

    std::vector<int> data = {1, 1};
    CHECK_HDF5(H5LTset_attribute_int(fileId, "/h5md", "version", data.data(), data.size()));

    CHECK_HDF5(H5LTset_attribute_string(fileId, "/h5md/author", "name", author.c_str()));

    std::string software = "espressopp";
    CHECK_HDF5(H5LTset_attribute_string(fileId, "/h5md/creator", "name", software.c_str()));

    std::string version = Version().info();
    CHECK_HDF5(H5LTset_attribute_string(fileId, "/h5md/creator", "version", version.c_str()));
}

void DumpH5MDParallel::writeBox(hid_t fileId)
{
    std::string groupName = "/particles/" + particleGroupName + "/box";
    auto group =
        CHECK_HDF5(H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    std::vector<int> dims = {3};
    CHECK_HDF5(
        H5LTset_attribute_int(fileId, groupName.c_str(), "dimension", dims.data(), dims.size()));

    auto boundaryType = H5Tcopy(H5T_C_S1);
    CHECK_HDF5(H5Tset_size(boundaryType, 8));
    CHECK_HDF5(H5Tset_strpad(boundaryType, H5T_STR_NULLPAD));
    std::vector<hsize_t> boundaryDims = {3};
    auto space = H5Screate_simple(int_c(boundaryDims.size()), boundaryDims.data(), nullptr);
    auto att = H5Acreate(group, "boundary", boundaryType, space, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_HDF5(H5Awrite(att, boundaryType, "periodicperiodicperiodic"));
    CHECK_HDF5(H5Aclose(att));
    CHECK_HDF5(H5Sclose(space));
    CHECK_HDF5(H5Tclose(boundaryType));

    std::vector<double> edges = {system_->bc->getBoxL()[0], system_->bc->getBoxL()[1],
                                 system_->bc->getBoxL()[2]};
    CHECK_HDF5(
        H5LTset_attribute_double(fileId, groupName.c_str(), "edges", edges.data(), edges.size()));

    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeId(hid_t fileId)
{
    using Datatype = int64_t;
    constexpr int64_t dimensions = 1;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + idDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->id());
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeType(hid_t fileId)
{
    using Datatype = int64_t;
    constexpr int64_t dimensions = 1;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + typeDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->type());
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeMass(hid_t fileId)
{
    using Datatype = double;
    constexpr int64_t dimensions = 1;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + massDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->mass());
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeQ(hid_t fileId)
{
    using Datatype = double;
    constexpr int64_t dimensions = 1;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + qDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->q());
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeGhost(hid_t fileId)
{
    using Datatype = int8_t;
    constexpr int64_t dimensions = 1;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + ghostDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->ghost());
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writePosition(hid_t fileId)
{
    using Datatype = double;
    constexpr int64_t dimensions = 3;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + positionDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->position()[0]);
        data.emplace_back(cit->position()[1]);
        data.emplace_back(cit->position()[2]);
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeVelocity(hid_t fileId)
{
    using Datatype = double;
    constexpr int64_t dimensions = 3;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + velocityDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->velocity()[0]);
        data.emplace_back(cit->velocity()[1]);
        data.emplace_back(cit->velocity()[2]);
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::writeForce(hid_t fileId)
{
    using Datatype = double;
    constexpr int64_t dimensions = 3;  ///< dimensions of the property

    std::string groupName = "/particles/" + particleGroupName + "/" + forceDataset;
    auto group = H5Gcreate(fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<Datatype> data;
    data.reserve(numLocalParticles * dimensions);
    for (iterator::CellListIterator cit(system_->storage->getRealCells()); !cit.isDone(); ++cit)
    {
        data.emplace_back(cit->force()[0]);
        data.emplace_back(cit->force()[1]);
        data.emplace_back(cit->force()[2]);
    }
    CHECK_EQUAL(int64_c(data.size()), numLocalParticles * dimensions);

    std::vector<hsize_t> localDims = {1, uint64_c(numLocalParticles), dimensions};
    std::vector<hsize_t> globalDims = {1, uint64_c(numTotalParticles), dimensions};

    std::string dataset_name = groupName + "/value";
    writeParallel(fileId, dataset_name, globalDims, localDims, data);

    std::vector<hsize_t> dims = {1};
    std::vector<int64_t> step = {0};
    std::vector<double> time = {0};
    std::string stepDataset = groupName + "/step";
    CHECK_HDF5(H5LTmake_dataset(fileId, stepDataset.c_str(), 1, dims.data(), typeToHDF5<int64_t>(),
                                step.data()));
    std::string timeDataset = groupName + "/time";
    CHECK_HDF5(H5LTmake_dataset(fileId, timeDataset.c_str(), 1, dims.data(), typeToHDF5<double>(),
                                time.data()));
    CHECK_HDF5(H5Gclose(group));
}

void DumpH5MDParallel::updateCache()
{
    boost::mpi::communicator world;
    comm = world;
    numProcesses = world.size();
    rank = world.rank();

    numLocalParticles = system_->storage->getNRealParticles();
    MPI_Allreduce(reinterpret_cast<const void*>(&numLocalParticles),
                  reinterpret_cast<void*>(&numTotalParticles), 1, MPI_INT64_T, MPI_SUM, comm);

    MPI_Exscan(&numLocalParticles, &particleOffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (rank == 0) particleOffset = 0;
}

void DumpH5MDParallel::dump()
{
    updateCache();

    auto info = MPI_INFO_NULL;

    auto plist = CHECK_HDF5(H5Pcreate(H5P_FILE_ACCESS));
    CHECK_HDF5(H5Pset_fapl_mpio(plist, comm, info));

    auto file_id = CHECK_HDF5(H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist));

    auto group1 =
        CHECK_HDF5(H5Gcreate(file_id, "/particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    std::string particleGroup = "/particles/" + particleGroupName;
    auto group2 = CHECK_HDF5(
        H5Gcreate(file_id, particleGroup.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    writeHeader(file_id);
    writeBox(file_id);
    if (dumpId) writeId(file_id);
    if (dumpType) writeType(file_id);
    if (dumpMass) writeMass(file_id);
    if (dumpQ) writeQ(file_id);
    if (dumpGhost) writeGhost(file_id);
    if (dumpPosition) writePosition(file_id);
    if (dumpVelocity) writeVelocity(file_id);
    if (dumpForce) writeForce(file_id);

    CHECK_HDF5(H5Gclose(group1));
    CHECK_HDF5(H5Gclose(group2));

    CHECK_HDF5(H5Fclose(file_id));
}

void DumpH5MDParallel::registerPython()
{
    using namespace espressopp::python;

    class_<DumpH5MDParallel>("io_DumpH5MDParallel", init<shared_ptr<System>, std::string>())
        .def_readwrite("author", &DumpH5MDParallel::author)
        .def_readwrite("particleGroupName", &DumpH5MDParallel::particleGroupName)
        .def_readonly("dumpId", &DumpH5MDParallel::dumpId)
        .def_readonly("dumpType", &DumpH5MDParallel::dumpType)
        .def_readonly("dumpMass", &DumpH5MDParallel::dumpMass)
        .def_readonly("dumpQ", &DumpH5MDParallel::dumpQ)
        .def_readonly("dumpGhost", &DumpH5MDParallel::dumpGhost)
        .def_readonly("dumpPosition", &DumpH5MDParallel::dumpPosition)
        .def_readonly("dumpVelocity", &DumpH5MDParallel::dumpVelocity)
        .def_readonly("dumpForce", &DumpH5MDParallel::dumpForce)
        .def_readwrite("idDataset", &DumpH5MDParallel::idDataset)
        .def_readwrite("typeDataset", &DumpH5MDParallel::typeDataset)
        .def_readwrite("massDataset", &DumpH5MDParallel::massDataset)
        .def_readwrite("qDataset", &DumpH5MDParallel::qDataset)
        .def_readwrite("ghostDataset", &DumpH5MDParallel::ghostDataset)
        .def_readwrite("positionDataset", &DumpH5MDParallel::positionDataset)
        .def_readwrite("velocityDataset", &DumpH5MDParallel::velocityDataset)
        .def_readwrite("forceDataset", &DumpH5MDParallel::forceDataset)
        .def("dump", &DumpH5MDParallel::dump);
}
}  // namespace io
}  // namespace espressopp