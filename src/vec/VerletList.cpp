/*
  Copyright (C) 2019-2021
      Max Planck Institute for Polymer Research & JGU Mainz
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

#include "Vectorization.hpp"
#include "VerletList.hpp"
#include "vec/storage/StorageVec.hpp"

#include "python.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

#include <atomic>

namespace espressopp
{
namespace vec
{
using namespace espressopp::iterator;

template <bool PACK_NEIGHBORS>
void rebuild_p_nc_pack_stencil(real const cutsq,
                               CellNeighborList const& cellNborList,
                               ParticleArray const& particleArray,
                               VerletList::NeighborList& neighborList);

LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/////////////////////////////////////////////////////////////////////////////////////////////////

/// cut is a cutoff (without skin)
VerletList::VerletList(std::shared_ptr<System> system, real _cut, bool rebuildVL)
    : SystemAccess(system)
{
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << _cut);

    if (!getSystem()->vectorization)
    {
        throw std::runtime_error("system has no vectorization");
    }
    vectorization = getSystem()->vectorization;

    if (!getSystem()->storage)
    {
        throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + getSystem()->getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    resetTimers();
    if (rebuildVL) rebuild();  // not called if exclutions are provided

    connect();
}

real VerletList::getVerletCutoff() { return cutVerlet; }

void VerletList::connect()
{
    // make a connection to vectorization to invoke rebuild on loadCells
    connectionResort =
        getSystem()->storage->onParticlesChanged.connect(std::bind(&VerletList::rebuild, this));
}

void VerletList::disconnect()
{
    // disconnect from System to avoid rebuild on resort
    connectionResort.disconnect();
}

void VerletList::rebuild()
{
    timer.reset();
    real currTime = timer.getElapsedTime();

    cutVerlet = cut + getSystem()->getSkin();
    cutsq = cutVerlet * cutVerlet;

    {
        neighborList.reset();
        const auto& cellnbrs = vectorization->neighborList;
        const auto& particles = vectorization->particles;

        if (particles.size())
        {
            rebuild_p_nc_pack_stencil<1>(cutsq, cellnbrs, particles, neighborList);
        }
    }
    timeRebuild += timer.getElapsedTime() - currTime;
    builds++;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

template <bool PACK_NEIGHBORS>
void rebuild_p_nc_pack_stencil(real const cutsq,
                               CellNeighborList const& cellNborList,
                               ParticleArray const& particleArray,
                               VerletList::NeighborList& neighborList)
{
    {
        const size_t* __restrict cellRange = particleArray.cellRange().data();
        const size_t* __restrict sizes = particleArray.sizes().data();

        // number of cells with neighbors
        const size_t numRealCells = cellNborList.numCells();
        const size_t numCells = particleArray.sizes().size();

        if (PACK_NEIGHBORS)
        {
            // real currTime = timer.getElapsedTime();

            size_t max_reserve = 0;
            for (size_t irow = 0; irow < numRealCells; irow++)
            {
                size_t c_reserve = 0;
                const size_t cell_nnbrs = cellNborList.numNeighbors(irow);
                for (size_t inbr = 0; inbr < cell_nnbrs; inbr++)
                {
                    const size_t cell_id = cellNborList.at(irow, inbr);
                    c_reserve += sizes[cell_id];
                }
                c_reserve = ESPP_FIT_TO_VECTOR_WIDTH(c_reserve);
                max_reserve = std::max(max_reserve, c_reserve);
            }

            if (max_reserve > neighborList.c_j.size())
            {
                neighborList.c_j.resize(max_reserve);
                neighborList.c_x.resize(max_reserve);
                neighborList.c_y.resize(max_reserve);
                neighborList.c_z.resize(max_reserve);
            }

            // timePack += timer.getElapsedTime() - currTime;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////

        const auto* __restrict pa_p_x = particleArray.p_x.data();
        const auto* __restrict pa_p_y = particleArray.p_y.data();
        const auto* __restrict pa_p_z = particleArray.p_z.data();
        const auto* __restrict pa_p_type = particleArray.type.data();

        auto* __restrict c_x_ptr = neighborList.c_x.data();
        auto* __restrict c_y_ptr = neighborList.c_y.data();
        auto* __restrict c_z_ptr = neighborList.c_z.data();
        auto* __restrict c_j_ptr = neighborList.c_j.data();

        int c_np_start = 0;
        for (size_t irow = 0; irow < numRealCells; irow++)
        {
            const size_t cell_id = cellNborList.cellId(irow);
            const size_t cell_nnbrs = cellNborList.numNeighbors(irow);
            const size_t cell_start = cellRange[cell_id];
            const size_t cell_data_end = cellRange[cell_id + 1];
            const size_t cell_size = sizes[cell_id];
            const size_t cell_end = cell_start + cell_size;

            int c_j_size = 0;
            if (PACK_NEIGHBORS)
            {
                int& nc_ctr = c_j_size;
                for (int inbr = 0; inbr < cell_nnbrs; inbr++)
                {
                    const int ncell_id = cellNborList.at(irow, inbr);
                    const int ncell_start = cellRange[ncell_id];
                    const int ncell_size = sizes[ncell_id];

                    int* __restrict c_j_ctr = c_j_ptr + nc_ctr;

#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
                    for (int ll = 0; ll < ncell_size; ll++)
                    {
                        c_j_ctr[ll] = ncell_start + ll;
                    }
                    nc_ctr += ncell_size;
                }

                /// padding
                {
                    const int last_nbr = cellNborList.at(irow, cell_nnbrs - 1);
                    const int padding = cellRange[last_nbr] + sizes[last_nbr];
                    const int pad_end = ESPP_FIT_TO_VECTOR_WIDTH(nc_ctr);
                    const int num_pad = pad_end - nc_ctr;
                    int* __restrict c_j_ctr = c_j_ptr + nc_ctr;

#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
                    for (size_t ll = 0; ll < num_pad; ll++)
                    {
                        c_j_ctr[ll] = padding;
                    }
                    nc_ctr += num_pad;
                }
                const int& c_end = nc_ctr;

                /// fill values
#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
                for (int ii = 0; ii < c_end; ii++)
                {
                    int p = c_j_ptr[ii];
                    c_x_ptr[ii] = pa_p_x[p];
                    c_y_ptr[ii] = pa_p_y[p];
                    c_z_ptr[ii] = pa_p_z[p];
                }
            }
            else
            {
                /// just count the total number of neighbors
                for (size_t inbr = 0; inbr < cell_nnbrs; inbr++)
                {
                    auto ncell_id = cellNborList.at(irow, inbr);
                    auto ncell_start = cellRange[ncell_id];
                    auto ncell_end = ncell_start + sizes[ncell_id];
                    c_j_size += ncell_end - ncell_start;
                }
            }

            /// Check whether there is enough space for all possible entries for
            /// this cell otherwise resize
            {
                const size_t c_np_size_max =
                    ESPP_FIT_TO_VECTOR_WIDTH((cell_size * (cell_size + c_j_size)));
                const size_t np_size_max = c_np_start + c_np_size_max;

                if (np_size_max > neighborList.nplist.size())
                {
                    neighborList.nplist.resize(2 * np_size_max);
                }
            }
            auto* __restrict nplist = neighborList.nplist.data();

            /// track counts of p and np for this cell
            int c_nplist_size = 0;

            for (size_t p = cell_start; p < cell_end; p++)
            {
                real p_x, p_y, p_z;
                p_x = pa_p_x[p];
                p_y = pa_p_y[p];
                p_z = pa_p_z[p];
                const int prev_c_nplist_size = c_nplist_size;

                // self-loop
                {
                    size_t ncell_id = cell_id;

                    size_t ncell_start = cellRange[ncell_id];
                    size_t ncell_data_end = cellRange[ncell_id + 1];
                    size_t ncell_end = ncell_start + sizes[ncell_id];
                    int* __restrict npptr = &(neighborList.nplist[c_np_start + c_nplist_size]);

                    {
                        int ll = 0;

#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
                        for (size_t np = ncell_start; np < ncell_data_end; np++)
                        {
                            {
                                real dist_x, dist_y, dist_z;
                                dist_x = p_x - pa_p_x[np];
                                dist_y = p_y - pa_p_y[np];
                                dist_z = p_z - pa_p_z[np];
                                const real distSqr =
                                    dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;

                                if (p < np && distSqr <= cutsq)
                                {
                                    npptr[ll++] = np;
                                }
                            }
                        }
                        c_nplist_size += ll;
                    }
                }

                size_t last_ncell = cell_id;
                size_t prev_nplist_size_nbrloop = c_nplist_size;

                /// neighbor-loop
                if (PACK_NEIGHBORS)
                {
                    auto* __restrict c_x_ptr = neighborList.c_x.data();
                    auto* __restrict c_y_ptr = neighborList.c_y.data();
                    auto* __restrict c_z_ptr = neighborList.c_z.data();
                    auto* __restrict c_j_ptr = neighborList.c_j.data();

                    int* __restrict npptr = &(neighborList.nplist[c_np_start + c_nplist_size]);
                    int ll = 0;

#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
                    for (int ii = 0; ii < c_j_size; ii++)
                    {
                        const real dist_x = p_x - c_x_ptr[ii];
                        const real dist_y = p_y - c_y_ptr[ii];
                        const real dist_z = p_z - c_z_ptr[ii];

                        const real distSqr = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;

                        int include_np = 1;
                        if (distSqr > cutsq) include_np = 0;
                        if (include_np) npptr[ll++] = c_j_ptr[ii];
                    }
                    c_nplist_size += ll;

                    // NOTE: last_ncell is not filled out correctly
                    // track the last cell to have a neighbor of p
                }
                else
                    for (size_t inbr = 0; inbr < cell_nnbrs; inbr++)
                    {
                        size_t ncell_id = cellNborList.at(irow, inbr);

                        size_t ncell_start = cellRange[ncell_id];
                        size_t ncell_data_end = cellRange[ncell_id + 1];
                        size_t ncell_end = ncell_start + sizes[ncell_id];
                        int* __restrict npptr = &(neighborList.nplist[c_np_start + c_nplist_size]);

                        {
                            int ll = 0;

#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
                            for (size_t np = ncell_start; np < ncell_data_end; np++)
                            {
                                {
                                    real dist_x, dist_y, dist_z;
                                    dist_x = p_x - pa_p_x[np];
                                    dist_y = p_y - pa_p_y[np];
                                    dist_z = p_z - pa_p_z[np];
                                    const real distSqr =
                                        dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;

                                    if (distSqr <= cutsq)
                                    {
                                        npptr[ll++] = np;
                                    }
                                }
                            }

                            c_nplist_size += ll;
                        }

                        // track the last cell to have a neighbor of p
                        if (c_nplist_size - prev_nplist_size_nbrloop) last_ncell = ncell_id;
                        prev_nplist_size_nbrloop = c_nplist_size;
                    }

                const int new_pairs = c_nplist_size - prev_c_nplist_size;

                if (new_pairs)
                {
                    // pad remaining part of list with stray neighbor particle
                    size_t num_rem = new_pairs % ESPP_VECTOR_WIDTH;
                    size_t num_pad = (num_rem > 0) * (ESPP_VECTOR_WIDTH - num_rem);
                    size_t padding = cellRange[last_ncell] + sizes[last_ncell];
                    for (size_t pad = 0; pad < num_pad; pad++)
                        nplist[c_np_start + (c_nplist_size++)] = padding;
                    neighborList.plist.push_back(p);

                    const int prange_start = c_np_start + prev_c_nplist_size;
                    const int prange_end = c_np_start + c_nplist_size;
                    neighborList.prange.push_back({prange_start, prange_end});
                }
            }
            c_np_start += c_nplist_size;
            // VEC_ASSERT_LEQ(c_np_start, neighborList.nplist.size())

            // TODO: require reduction for num_pairs or atomic update
            // num_pairs += c_nplist_size;

            size_t max_type_cell = 0;
            for (size_t p = cell_start; p < cell_end; p++)
                max_type_cell = std::max(max_type_cell, pa_p_type[p]);
            neighborList.max_type = std::max(neighborList.max_type, max_type_cell);
        }
        neighborList.num_pairs = c_np_start;

#if 0
      {
        // real currTime = timer.getElapsedTime();
#if 1
        {
          for(size_t irow=0; irow < numRealCells; irow++) f_cell(irow);
        }
#else
        {
          // executeInScheduler(f_cell);
          using hpx::parallel::execution::par;
          using hpx::parallel::execution::seq;
          using hpx::parallel::for_loop;
          utils::parallelForLoop(size_t(0), numRealCells, f_cell);
        }
#endif
        // timeExecute += timer.getElapsedTime() - currTime;
      }
#endif

#if 0
      {
        int num_pairs_check = 0;
        for(size_t irow=0; irow < numRealCells; irow++)
        {
          const auto& crange = neighborList.crange[irow];
          for(int ip=crange.first; ip<crange.second; ip++){
            const auto& prange = neighborList.prange[ip];
            num_pairs_check += (prange.second - prange.first);
          }
        }
        if(num_pairs_check!=num_pairs){
          VEC_THROW_EXCEPTION(hpx::assertion_failure,"VerletList::rebuild_p_nc_pack_stencil",
            "num_pairs_check: Expected "<<num_pairs<<" got "<<num_pairs_check);
        }
      }
#endif

#if 0
      {
        const auto* __restrict pa_p_type = particleArray.type.data();

        size_t max_type_check = 0;
        for(size_t irow=0; irow < numRealCells; irow++)
        {
          const size_t cell_id       = cellNborList.cellId(irow);
          const size_t cell_start    = cellRange[cell_id];
          const size_t cell_size     = sizes[cell_id];
          const size_t cell_end      = cell_start + cell_size;
          for(size_t p=cell_start; p < cell_end; p++)
            max_type_check = std::max(max_type_check, pa_p_type[p]);
        }
        if(max_type_check!=max_type){
          VEC_THROW_EXCEPTION(hpx::assertion_failure,"VerletList::rebuild_p_nc_pack_stencil",
            "max_type_check: Expected "<<max_type<<" got "<<max_type_check);
        }
      }
#endif

#if 0
      {
        const size_t nplist_reserve = neighborList.nplist.size();
        if((num_pairs>nplist_reserve)&&(nplist_reserve>0)) {
          LOG4ESPP_WARN(theLogger,"Reserve size exceeded. "
            "Expected "<< nplist_reserve<<". Got "<<num_pairs);
        }
        std::cout << "nplist_alloc: " << nplist_alloc << "num_pairs: " << num_pairs
          << " wasted: " << double(nplist_alloc-num_pairs)/double(nplist_alloc) << std::endl;
      }
#endif
    }

    // LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
    //              << " local size = " << vlPairs.size());
}

/////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
  void VerletList::checkPair(Particle& pt1, Particle& pt2)
  {

    Real3D d = pt1.position() - pt2.position();
    real distsq = d.sqr();

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id()
                   << " @ " << pt1.position()
       << " - p2: " << pt2.id() << " @ " << pt2.position()
       << " -> distsq = " << distsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (both directions)
    if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;

    max_type = std::max(max_type, (std::max(pt1.type(), pt2.type())) );
    vlPairs.add(pt1, pt2); // add pair to Verlet List
  }
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////

int VerletList::totalSize() const
{
    System& system = getSystemRef();
    int size = localSize();
    int allsize;

    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
}

int VerletList::localSize() const { return neighborList.num_pairs; }

bool VerletList::exclude(longint pid1, longint pid2)
{
    throw std::runtime_error("Exclusions in VerletList not implemented.");
    // exList.insert(std::make_pair(pid1, pid2));
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

VerletList::~VerletList()
{
    LOG4ESPP_INFO(theLogger, "~VerletList");

    if (!connectionResort.connected())
    {
        connectionResort.disconnect();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void VerletList::resetTimers() { timeRebuild = 0.0; }

void VerletList::loadTimers(real* t)
{
    t[0] = timeRebuild;
    ;
}

static boost::python::object wrapGetTimers(class VerletList* obj)
{
    real tms[1];
    obj->loadTimers(tms);
    return boost::python::make_tuple(tms[0]);
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void VerletList::registerPython()
{
    using namespace espressopp::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2) = &VerletList::exclude;

    class_<VerletList, std::shared_ptr<VerletList> >("vec_VerletList",
                                                     init<std::shared_ptr<System>, real, bool>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
        .def("totalSize", &VerletList::totalSize)
        .def("localSize", &VerletList::localSize)
        // .def("getPair", &VerletList::getPair)
        // .def("exclude", pyExclude)
        .def("rebuild", &VerletList::rebuild)
        // .def("rebuildPairs", &VerletList::rebuildPairs)
        .def("connect", &VerletList::connect)
        .def("disconnect", &VerletList::disconnect)
        .def("getVerletCutoff", &VerletList::getVerletCutoff)
        .def("resetTimers", &VerletList::resetTimers)
        .def("getTimers", wrapGetTimers);
}

}  // namespace vec
}  // namespace espressopp
