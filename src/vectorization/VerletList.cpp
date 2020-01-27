/*
  Copyright (C) 2019
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

#include "python.hpp"
#include "VerletList.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "Vectorization.hpp"

namespace espressopp { namespace vectorization {

  using namespace espressopp::iterator;

  LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

  // cut is a cutoff (without skin)
  VerletList::VerletList(shared_ptr<System> system, shared_ptr<Vectorization> vec,
    real _cut, bool rebuildVL, int build_order)
    : SystemAccess(system), vec(vec), build_order(build_order)
  {
    std::string build_order_str[] = {"c->p->nc->np","c->nc->p->np"};
    LOG4ESPP_INFO(theLogger,"Using vectorized verlet list with build order " << build_order_str[build_order]);
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << _cut);

    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;
    max_type = 0;

    resetTimers();
    if (rebuildVL) rebuild(); // not called if exclutions are provided

    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletList::rebuild, this));
  }

  real VerletList::getVerletCutoff(){
    return cutVerlet;
  }

  void VerletList::connect()
  {

  // make a connection to System to invoke rebuild on resort
  connectionResort = getSystem()->storage->onParticlesChanged.connect(
      boost::bind(&VerletList::rebuild, this));
  }

  void VerletList::disconnect()
  {

  // disconnect from System to avoid rebuild on resort
  connectionResort.disconnect();
  }

  /*-------------------------------------------------------------*/

  ParticleArray& VerletList::getParticleArray() {
    return vec->getParticleArray();
  }

  /*-------------------------------------------------------------*/

  void VerletList::rebuildPairs()
  {
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;

    vlPairs.clear();

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      checkPair(*it->first, *it->second);
      LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " and " << it->second->id());
    }
  }

  void VerletList::rebuild()
  {
    timer.reset();
    real currTime = timer.getElapsedTime();

    //real cutVerlet = cut + getSystem() -> getSkin();
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;

    vlPairs.clear();

    num_pairs = 0;
    neighborList.reset();

    const bool VEC_MODE_AOS = getParticleArray().mode_aos();
    const bool PACK_NEIGHBORS = true;

    if(build_order==0)
    {
      if(VEC_MODE_AOS)
      {
        if(PACK_NEIGHBORS)
          rebuild_p_nc_pack_stencil<true,true>();
        else
          rebuild_p_nc_pack_stencil<true,false>();
      }
      else
      {
        if(PACK_NEIGHBORS)
          rebuild_p_nc_pack_stencil<false,true>();
        else
          rebuild_p_nc_pack_stencil<false,false>();
      }
    }
    else if(build_order==1)
      rebuild_nc_p();
    else{
      LOG4ESPP_ERROR(theLogger,"Invalid build_order: "<<build_order);
      throw std::runtime_error("Invalid build_order");
    }

    timeRebuild += timer.getElapsedTime() - currTime;
    builds++;
  }

  /*-------------------------------------------------------------*/

  template< bool VEC_MODE_AOS >
  void VerletList::rebuild_p_nc()
  {
    {
      const auto& cellNborList            = vec->getNeighborList();
      const auto& particleArray           = vec->getParticleArray();
      const size_t* __restrict cellRange = &(particleArray.cellRange()[0]);
      const size_t* __restrict sizes      = &(particleArray.sizes()[0]);

      const auto* __restrict position  = &(particleArray.position[0]);
      const auto* __restrict pa_p_x    = &(particleArray.p_x[0]);
      const auto* __restrict pa_p_y    = &(particleArray.p_y[0]);
      const auto* __restrict pa_p_z    = &(particleArray.p_z[0]);
      const auto* __restrict pa_p_type = &(particleArray.type[0]);

      // number of cells with neighbors
      const size_t numRealCells = cellNborList.numCells();
      const size_t numCells = particleArray.sizes().size();

      neighborList.clist.resize(numRealCells);
      neighborList.crange.resize(numRealCells);

      //// preallocate space
      size_t plist_reserve = 0;
      for(size_t irow=0; irow<numRealCells; irow++)
      {
        plist_reserve += sizes[cellNborList.cellId(irow)];
      }
      neighborList.plist.reserve(plist_reserve);
      neighborList.prange.reserve(plist_reserve);

      size_t max_nneighbors = cellNborList.maxNumNeighbors();
      size_t max_cell_size = 0, nplist_reserve;
      for(size_t icell=0; icell<numCells; icell++)
        max_cell_size = std::max(max_cell_size, sizes[icell]);
      size_t max_pairs_per_cell = max_cell_size*max_cell_size*(max_nneighbors+1);
      nplist_reserve = neighborList.nplist.size(); // re-use previous allocation

      size_t prev_num_pairs = num_pairs;
      for(size_t irow=0; irow < numRealCells; irow++)
      {
        ///////////////////////////////////////////////////////////////////////////////////////////
        /// Check whether there is enough space for next cells
        if(num_pairs + max_pairs_per_cell > (neighborList.nplist.size()))
        {
          Int3D cellgrid = getSystem()->storage->getInt3DCellGrid();
          const size_t threshold = std::max(cellgrid[0],std::max(cellgrid[1],cellgrid[2]));
          size_t cells_left = numRealCells-irow, add_size;
          if(cells_left <= threshold)
            // over-estimate with few cells left
            add_size = max_pairs_per_cell*(numRealCells-irow);
          else
            // under-estimate with many cells left
            add_size = max_pairs_per_cell*(numRealCells-irow) / threshold;
          size_t old_size = neighborList.nplist.size();
          size_t new_size = ESPP_FIT_TO_VECTOR_WIDTH(old_size+add_size);
          neighborList.nplist.resize(new_size);
          LOG4ESPP_INFO(theLogger,"Resizing neighbor list. num_pairs="<<num_pairs
            <<" old_size="<<old_size<<" new_size="<<new_size);
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        size_t  cell_id         = cellNborList.cellId(irow);
        size_t  cell_nnbrs      = cellNborList.numNeighbors(irow);
        size_t  cell_start      = cellRange[cell_id];
        size_t  cell_data_end   = cellRange[cell_id+1];
        size_t  cell_size       = sizes[cell_id];
        size_t  cell_end        = cell_start + cell_size;

        for(size_t p=cell_start; p < cell_end; p++)
        {
          real p_x, p_y, p_z;
          if(VEC_MODE_AOS)
          {
            p_x = position[p].x;
            p_y = position[p].y;
            p_z = position[p].z;
          }
          else
          {
            p_x = pa_p_x[p];
            p_y = pa_p_y[p];
            p_z = pa_p_z[p];
          }

          // self-loop
          {
            size_t ncell_id = cell_id;

            size_t  ncell_start     = cellRange[ncell_id];
            size_t  ncell_data_end  = cellRange[ncell_id+1];
            size_t  ncell_end       = ncell_start + sizes[ncell_id];
            int* __restrict npptr =  &(neighborList.nplist[num_pairs]);

            {
              int ll=0;

              #pragma vector always
              #pragma vector aligned
              #pragma ivdep
              for(size_t np = ncell_start; np<ncell_data_end; np++)
              {
                {
                  int include_np = 1;
                  real dist_x, dist_y, dist_z;
                  if(VEC_MODE_AOS)
                  {
                    dist_x   = p_x - position[np].x;
                    dist_y   = p_y - position[np].y;
                    dist_z   = p_z - position[np].z;
                  }
                  else
                  {
                    dist_x   = p_x - pa_p_x[np];
                    dist_y   = p_y - pa_p_y[np];
                    dist_z   = p_z - pa_p_z[np];
                  }

                  const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                  if(p>=np) include_np = 0;

                  if(distSqr > cutsq) include_np = 0;
                  if(include_np) npptr[ll++] = np;
                }
              }
              num_pairs += ll;
            }
          }

          // neighbor-loop
          size_t last_ncell = cell_id, prev_num_pairs_2 = num_pairs;
          for(size_t inbr=0; inbr<cell_nnbrs; inbr++)
          {
            size_t ncell_id = cellNborList.at(irow,inbr);

            size_t  ncell_start     = cellRange[ncell_id];
            size_t  ncell_data_end  = cellRange[ncell_id+1];
            size_t  ncell_end       = ncell_start + sizes[ncell_id];
            int* __restrict npptr =  &(neighborList.nplist[num_pairs]);

            {
              int ll=0;

              #pragma vector always
              #pragma vector aligned
              #pragma ivdep
              for(size_t np = ncell_start; np<ncell_data_end; np++)
              {
                {
                  int include_np = 1;
                  real dist_x, dist_y, dist_z;
                  if(VEC_MODE_AOS)
                  {
                    dist_x   = p_x - position[np].x;
                    dist_y   = p_y - position[np].y;
                    dist_z   = p_z - position[np].z;
                  }
                  else
                  {
                    dist_x   = p_x - pa_p_x[np];
                    dist_y   = p_y - pa_p_y[np];
                    dist_z   = p_z - pa_p_z[np];
                  }

                  const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                  if(distSqr > cutsq) include_np = 0;
                  if(include_np) npptr[ll++] = np;
                }
              }

              num_pairs += ll;
            }

            // track the last cell to have a neighbor of p
            if(num_pairs - prev_num_pairs_2)
              last_ncell = ncell_id;
            prev_num_pairs_2 = num_pairs;
          }

          size_t new_pairs = (num_pairs - prev_num_pairs);

          if(new_pairs)
          {
            // pad remaining part of list with stray neighbor particle
            size_t num_rem = num_pairs % ESPP_VECTOR_WIDTH;
            size_t num_pad = (num_rem > 0) * (ESPP_VECTOR_WIDTH - num_rem);
            size_t padding = cellRange[last_ncell]+sizes[last_ncell];
            for(size_t pad = 0; pad < num_pad; pad++)
              neighborList.nplist[num_pairs++] = padding;
            neighborList.plist.push_back(p);
            neighborList.prange.push_back(num_pairs);
          }
          prev_num_pairs = num_pairs;
        }

        if(VEC_MODE_AOS)
          for(size_t p=cell_start; p < cell_end; p++)
            max_type = std::max(max_type, position[p].t);
        else
          for(size_t p=cell_start; p < cell_end; p++)
            max_type = std::max(max_type, pa_p_type[p]);
      }

      if(num_pairs>nplist_reserve) {
        LOG4ESPP_WARN(theLogger,"Reserve size exceeded. "
          "Expected "<< nplist_reserve<<". Got "<<num_pairs);
      }
    }

    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlPairs.size());
  }

  /*-------------------------------------------------------------*/

  template< bool VEC_MODE_AOS, bool PACK_NEIGHBORS >
  void VerletList::rebuild_p_nc_pack_stencil()
  {
    {
      const auto& cellNborList            = vec->getNeighborList();
      const auto& particleArray           = vec->getParticleArray();
      const size_t* __restrict cellRange = &(particleArray.cellRange()[0]);
      const size_t* __restrict sizes      = &(particleArray.sizes()[0]);

      const auto* __restrict position  = &(particleArray.position[0]);
      const auto* __restrict pa_p_x    = &(particleArray.p_x[0]);
      const auto* __restrict pa_p_y    = &(particleArray.p_y[0]);
      const auto* __restrict pa_p_z    = &(particleArray.p_z[0]);
      const auto* __restrict pa_p_type = &(particleArray.type[0]);

      // number of cells with neighbors
      const size_t numRealCells = cellNborList.numCells();
      const size_t numCells = particleArray.sizes().size();

      neighborList.clist.resize(numRealCells);
      neighborList.crange.resize(numRealCells);

      //// preallocate space
      size_t plist_reserve = 0;
      for(size_t irow=0; irow<numRealCells; irow++)
      {
        plist_reserve += sizes[cellNborList.cellId(irow)];
      }
      neighborList.plist.reserve(plist_reserve);
      neighborList.prange.reserve(plist_reserve);

      size_t max_nneighbors = cellNborList.maxNumNeighbors();
      size_t max_cell_size = 0, nplist_reserve;
      for(size_t icell=0; icell<numCells; icell++)
        max_cell_size = std::max(max_cell_size, sizes[icell]);
      size_t max_pairs_per_cell = max_cell_size*max_cell_size*(max_nneighbors+1);
      nplist_reserve = neighborList.nplist.size(); // re-use previous allocation

      /////////////////////////////////////////////////////////////////////////////////////////////
      ///// pack neighbors together

      std::vector<int> c_range;

      if(PACK_NEIGHBORS)
      {
        c_range.reserve(numRealCells+1);
        size_t vec_max_cell_size = ((max_cell_size+ESPP_VECTOR_WIDTH-1)/ESPP_VECTOR_WIDTH)*ESPP_VECTOR_WIDTH;
        size_t c_reserve = numRealCells*max_nneighbors*vec_max_cell_size;
        if(c_reserve > c_j.size())
        {
          c_j.resize(c_reserve);
          c_x.resize(c_reserve);
          c_y.resize(c_reserve);
          c_z.resize(c_reserve);
        }
      }

      auto* __restrict c_x_ptr = &(c_x[0]);
      auto* __restrict c_y_ptr = &(c_y[0]);
      auto* __restrict c_z_ptr = &(c_z[0]);
      auto* __restrict c_j_ptr = &(c_j[0]);
      int c_j_max = 0;

      if(PACK_NEIGHBORS)
      {
        c_range.push_back(c_j_max);
        for(size_t irow=0; irow<numRealCells; irow++)
        {
          const size_t cell_nnbrs = cellNborList.numNeighbors(irow);
          for(size_t inbr=0; inbr<cell_nnbrs; inbr++)
          {
            size_t cell_id       = cellNborList.at(irow,inbr);
            size_t cell_start    = cellRange[cell_id];
            size_t cell_size     = sizes[cell_id];
            size_t cell_end      = cell_start + cell_size;
            int* __restrict c_j_ctr = c_j_ptr + c_j_max;
            int ll = 0;

            #pragma vector always
            #pragma vector aligned
            #pragma ivdep
            for(size_t ll=0; ll<cell_size; ll++)
            {
              c_j_ctr[ll] = cell_start+ll;
            }
            c_j_max += cell_size;
          }

          // padding
          {
            int last_nbr = cellNborList.at(irow,cell_nnbrs-1);
            int padding = cellRange[last_nbr]+sizes[last_nbr];
            int pad_end = ((c_j_max+ESPP_VECTOR_WIDTH-1)/ESPP_VECTOR_WIDTH)*ESPP_VECTOR_WIDTH;
            int num_padding  = pad_end - c_j_max;
            int* __restrict c_j_ctr = c_j_ptr + c_j_max;

            #pragma vector always
            #pragma vector aligned
            #pragma ivdep
            for(size_t ll=0; ll<num_padding; ll++)
            {
              c_j_ctr[ll] = padding;
            }
            c_j_max += num_padding;
          }
          c_range.push_back(c_j_max);
        }
        if(c_j_max>c_j.size()) throw std::runtime_error("rebuild_p_nc_pack_stencil: Reserve size exceeded.");

        /// fill values
        #pragma vector always
        #pragma vector aligned
        #pragma ivdep
        for(int ii=0; ii<c_j_max; ii++)
        {
          int p = c_j_ptr[ii];
          if(VEC_MODE_AOS)
          {
            c_x_ptr[ii] = position[p].x;
            c_y_ptr[ii] = position[p].y;
            c_z_ptr[ii] = position[p].z;
          }
          else
          {
            c_x_ptr[ii] = pa_p_x[p];
            c_y_ptr[ii] = pa_p_y[p];
            c_z_ptr[ii] = pa_p_z[p];
          }
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////

      size_t prev_num_pairs = num_pairs;
      for(size_t irow=0; irow < numRealCells; irow++)
      {
        ///////////////////////////////////////////////////////////////////////////////////////////
        /// Check whether there is enough space for next cells
        if(num_pairs + max_pairs_per_cell > (neighborList.nplist.size()))
        {
          Int3D cellgrid = getSystem()->storage->getInt3DCellGrid();
          const size_t threshold = std::max(cellgrid[0],std::max(cellgrid[1],cellgrid[2]));
          size_t cells_left = numRealCells-irow, add_size;
          if(cells_left <= threshold)
            // over-estimate with few cells left
            add_size = max_pairs_per_cell*(numRealCells-irow);
          else
            // under-estimate with many cells left
            add_size = max_pairs_per_cell*(numRealCells-irow) / threshold;
          size_t old_size = neighborList.nplist.size();
          size_t new_size = ESPP_FIT_TO_VECTOR_WIDTH(old_size+add_size);
          neighborList.nplist.resize(new_size);
          LOG4ESPP_INFO(theLogger,"Resizing neighbor list. num_pairs="<<num_pairs
            <<" old_size="<<old_size<<" new_size="<<new_size);
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        size_t  cell_id       = cellNborList.cellId(irow);
        size_t  cell_nnbrs    = cellNborList.numNeighbors(irow);
        size_t  cell_start    = cellRange[cell_id];
        size_t  cell_data_end = cellRange[cell_id+1];
        size_t  cell_size     = sizes[cell_id];
        size_t  cell_end      = cell_start + cell_size;

        for(size_t p=cell_start; p < cell_end; p++)
        {
          real p_x, p_y, p_z;
          if(VEC_MODE_AOS)
          {
            p_x = position[p].x;
            p_y = position[p].y;
            p_z = position[p].z;
          }
          else
          {
            p_x = pa_p_x[p];
            p_y = pa_p_y[p];
            p_z = pa_p_z[p];
          }

          // self-loop
          {
            size_t ncell_id = cell_id;

            size_t  ncell_start    = cellRange[ncell_id];
            size_t  ncell_data_end = cellRange[ncell_id+1];
            size_t  ncell_end      = ncell_start + sizes[ncell_id];
            int* __restrict npptr  = &(neighborList.nplist[num_pairs]);

            {
              int ll=0;

              #pragma vector always
              #pragma vector aligned
              #pragma ivdep
              for(size_t np = ncell_start; np<ncell_data_end; np++)
              {
                {
                  int include_np = 1;
                  real dist_x, dist_y, dist_z;
                  if(VEC_MODE_AOS)
                  {
                    dist_x   = p_x - position[np].x;
                    dist_y   = p_y - position[np].y;
                    dist_z   = p_z - position[np].z;
                  }
                  else
                  {
                    dist_x   = p_x - pa_p_x[np];
                    dist_y   = p_y - pa_p_y[np];
                    dist_z   = p_z - pa_p_z[np];
                  }

                  const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                  if(p>=np) include_np = 0;

                  if(distSqr > cutsq) include_np = 0;
                  if(include_np) npptr[ll++] = np;
                }
              }
              num_pairs += ll;
            }
          }

          // neighbor-loop
          size_t last_ncell = cell_id, prev_num_pairs_2 = num_pairs;

          if(PACK_NEIGHBORS)
          {
            const int c_start = c_range[irow];
            const int c_end   = c_range[irow+1];
            int* __restrict npptr =  &(neighborList.nplist[num_pairs]);
            int ll=0;

            #pragma vector always
            #pragma vector aligned
            #pragma ivdep
            for(int ii=c_start; ii<c_end; ii++)
            {
              const real dist_x = p_x - c_x_ptr[ii];
              const real dist_y = p_y - c_y_ptr[ii];
              const real dist_z = p_z - c_z_ptr[ii];

              const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

              int include_np = 1;
              if(distSqr > cutsq) include_np = 0;
              if(include_np) npptr[ll++] = c_j_ptr[ii];
            }
            num_pairs += ll;
          }
          else
          for(size_t inbr=0; inbr<cell_nnbrs; inbr++)
          {
            size_t ncell_id = cellNborList.at(irow,inbr);

            size_t  ncell_start     = cellRange[ncell_id];
            size_t  ncell_data_end  = cellRange[ncell_id+1];
            size_t  ncell_end       = ncell_start + sizes[ncell_id];
            int* __restrict npptr =  &(neighborList.nplist[num_pairs]);

            {
              int ll=0;

              #pragma vector always
              #pragma vector aligned
              #pragma ivdep
              for(size_t np = ncell_start; np<ncell_data_end; np++)
              {
                {
                  int include_np = 1;
                  real dist_x, dist_y, dist_z;
                  if(VEC_MODE_AOS)
                  {
                    dist_x = p_x - position[np].x;
                    dist_y = p_y - position[np].y;
                    dist_z = p_z - position[np].z;
                  }
                  else
                  {
                    dist_x = p_x - pa_p_x[np];
                    dist_y = p_y - pa_p_y[np];
                    dist_z = p_z - pa_p_z[np];
                  }

                  const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                  if(distSqr > cutsq) include_np = 0;
                  if(include_np) npptr[ll++] = np;
                }
              }

              num_pairs += ll;
            }

            // track the last cell to have a neighbor of p
            if(num_pairs - prev_num_pairs_2)
              last_ncell = ncell_id;
            prev_num_pairs_2 = num_pairs;
          }

          size_t new_pairs = (num_pairs - prev_num_pairs);

          if(new_pairs)
          {
            // pad remaining part of list with stray neighbor particle
            size_t num_rem = num_pairs % ESPP_VECTOR_WIDTH;
            size_t num_pad = (num_rem > 0) * (ESPP_VECTOR_WIDTH - num_rem);
            size_t padding = cellRange[last_ncell]+sizes[last_ncell];
            for(size_t pad = 0; pad < num_pad; pad++)
              neighborList.nplist[num_pairs++] = padding;
            neighborList.plist.push_back(p);
            neighborList.prange.push_back(num_pairs);
          }
          prev_num_pairs = num_pairs;
        }

        if(VEC_MODE_AOS)
          for(size_t p=cell_start; p < cell_end; p++)
            max_type = std::max(max_type, position[p].t);
        else
          for(size_t p=cell_start; p < cell_end; p++)
            max_type = std::max(max_type, pa_p_type[p]);
      }

      if(num_pairs>nplist_reserve) {
        LOG4ESPP_WARN(theLogger,"Reserve size exceeded. "
          "Expected "<< nplist_reserve<<". Got "<<num_pairs);
      }
    }

    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlPairs.size());
  }

  /*-------------------------------------------------------------*/

  void VerletList::rebuild_nc_p()
  {
    throw std::runtime_error("VerletList::rebuild_nc_p(): Function not implemented.");
    {
      auto& cellNborList  = vec->getNeighborList();
      auto& particleArray = vec->getParticleArray();
      auto const& cellRange = particleArray.cellRange();
      auto const& sizes     = particleArray.sizes();

      #if defined(ESPP_AOS)
        const auto* __restrict position  = &(particleArray.position[0]);
      #else
        const auto* __restrict pa_p_x    = &(particleArray.p_x[0]);
        const auto* __restrict pa_p_y    = &(particleArray.p_y[0]);
        const auto* __restrict pa_p_z    = &(particleArray.p_z[0]);
        const auto* __restrict pa_p_type = &(particleArray.type[0]);
      #endif

      // number of cells with neighbors
      const size_t numRealCells = cellNborList.numCells();
      const size_t max_nneighbors = cellNborList.maxNumNeighbors();

      neighborList.clist.resize(numRealCells);
      neighborList.crange.resize(numRealCells);

      //// preallocate space
      size_t plist_reserve = 0;
      for(size_t irow=0; irow<numRealCells; irow++)
      {
        plist_reserve += sizes[cellNborList.cellId(irow)];
      }
      plist_reserve *= (max_nneighbors+1);
      neighborList.plist.reserve(plist_reserve);
      neighborList.prange.reserve(plist_reserve);

      size_t max_cell_size = 0, nplist_reserve;
      for(size_t icell=0; icell<sizes.size(); icell++)
        max_cell_size = std::max(max_cell_size, sizes[icell]);
      size_t max_pairs_per_cell = max_cell_size*max_cell_size*(max_nneighbors+1)*(max_nneighbors+1);
      nplist_reserve = neighborList.nplist.size(); // re-use previous allocation

      size_t prev_num_pairs = num_pairs;
      for(size_t irow=0; irow < numRealCells; irow++)
      {
        ///////////////////////////////////////////////////////////////////////////////////////////
        /// Check whether there is enough space for next cells
        if(num_pairs + max_pairs_per_cell > (neighborList.nplist.size()))
        {
          Int3D cellgrid = getSystem()->storage->getInt3DCellGrid();
          const size_t threshold = std::max(cellgrid[0],std::max(cellgrid[1],cellgrid[2]));
          size_t cells_left = numRealCells-irow, add_size;
          if(cells_left <= threshold)
            // over-estimate with few cells left
            add_size = max_pairs_per_cell*(numRealCells-irow);
          else
            // under-estimate with many cells left
            add_size = max_pairs_per_cell*(numRealCells-irow) / threshold;
          size_t old_size = neighborList.nplist.size();
          size_t new_size = ESPP_FIT_TO_VECTOR_WIDTH(old_size+add_size);
          neighborList.nplist.resize(new_size);
          LOG4ESPP_ERROR(theLogger,"Resizing neighbor list. num_pairs="<<num_pairs
            <<" old_size="<<old_size<<" new_size="<<new_size);
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        size_t  cell_id         = cellNborList.cellId(irow);
        size_t  cell_nnbrs      = cellNborList.numNeighbors(irow);
        size_t  cell_start      = cellRange[cell_id];
        size_t  cell_data_end   = cellRange[cell_id+1];
        size_t  cell_size       = sizes[cell_id];
        size_t  cell_end        = cell_start + cell_size;

        // self-loop
        {
          size_t ncell_id = cell_id;

          size_t  ncell_start     = cellRange[ncell_id];
          size_t  ncell_data_end  = cellRange[ncell_id+1];
          size_t  ncell_end       = ncell_start + particleArray.sizes()[ncell_id];

          for(size_t p=cell_start; p < cell_end; p++)
          {
            #if defined(ESPP_AOS)
              const real  p_x    = position[p].x;
              const real  p_y    = position[p].y;
              const real  p_z    = position[p].z;
            #else
              const real  p_x    = pa_p_x[p];
              const real  p_y    = pa_p_y[p];
              const real  p_z    = pa_p_z[p];
            #endif

            {
              int* __restrict npptr = &(neighborList.nplist[num_pairs]);
              int ll=0;

              #pragma vector always
              #pragma vector aligned
              #pragma ivdep
              for(size_t np = ncell_start; np<ncell_data_end; np++)
              {
                {
                  int include_np = 1;
                #if defined(ESPP_AOS)
                  const real dist_x   = p_x - position[np].x;
                  const real dist_y   = p_y - position[np].y;
                  const real dist_z   = p_z - position[np].z;
                #else
                  const real dist_x   = p_x - pa_p_x[np];
                  const real dist_y   = p_y - pa_p_y[np];
                  const real dist_z   = p_z - pa_p_z[np];
                #endif

                  const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                  if(p>=np) include_np = 0;

                  if(distSqr > cutsq) include_np = 0;
                  // if(np>=ncell_end) include_np = 0;
                  if(include_np) npptr[ll++] = np;
                }
              }
              num_pairs += ll;
              if(ll>0)
              {
                size_t num_rem = num_pairs % ESPP_VECTOR_WIDTH;
                size_t num_pad = (num_rem > 0) * (ESPP_VECTOR_WIDTH - num_rem);
                size_t padding = ncell_end;
                for(size_t pad = 0; pad < num_pad; pad++)
                  neighborList.nplist[num_pairs++] = padding;
                neighborList.plist.push_back(p);
                neighborList.prange.push_back(num_pairs);
              }
            }
          }
        }

        for(size_t inbr=0; inbr<cell_nnbrs; inbr++)
        {
          size_t ncell_id = cellNborList.at(irow,inbr);

          size_t  ncell_start     = cellRange[ncell_id];
          size_t  ncell_data_end  = cellRange[ncell_id+1];
          size_t  ncell_end       = ncell_start + particleArray.sizes()[ncell_id];

          for(size_t p=cell_start; p < cell_end; p++)
          {
            #if defined(ESPP_AOS)
              const real  p_x    = particleArray.position[p].x;
              const real  p_y    = particleArray.position[p].y;
              const real  p_z    = particleArray.position[p].z;
            #else
              const real  p_x    = pa_p_x[p];
              const real  p_y    = pa_p_y[p];
              const real  p_z    = pa_p_z[p];
            #endif

            {
              int* __restrict npptr = &(neighborList.nplist[num_pairs]);
              int ll=0;

              #pragma vector always
              #pragma vector aligned
              #pragma ivdep
              for(size_t np = ncell_start; np<ncell_data_end; np++)
              {
                {
                  int include_np = 1;
                #if defined(ESPP_AOS)
                  const real dist_x   = p_x - particleArray.position[np].x;
                  const real dist_y   = p_y - particleArray.position[np].y;
                  const real dist_z   = p_z - particleArray.position[np].z;
                #else
                  const real dist_x   = p_x - pa_p_x[np];
                  const real dist_y   = p_y - pa_p_y[np];
                  const real dist_z   = p_z - pa_p_z[np];
                #endif

                  const real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

                  if(distSqr > cutsq) include_np = 0;
                  if(include_np) npptr[ll++] = np;
                }
              }
              num_pairs += ll;
              if(ll>0)
              {
                size_t num_rem = num_pairs % ESPP_VECTOR_WIDTH;
                size_t num_pad = (num_rem > 0) * (ESPP_VECTOR_WIDTH - num_rem);
                size_t padding = ncell_end;
                for(size_t pad = 0; pad < num_pad; pad++)
                  neighborList.nplist[num_pairs++] = padding;
                neighborList.plist.push_back(p);
                neighborList.prange.push_back(num_pairs);
              }
            }
          }
        }

        for(size_t p=cell_start; p < cell_end; p++)
        {
         #if defined(ESPP_AOS)
          max_type = std::max(max_type, position[p].t);
         #else
          max_type = std::max(max_type, pa_p_type[p]);
         #endif
        }
      }
    }
  }

  /*-------------------------------------------------------------*/

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

    max_type = std::max(max_type, std::max(pt1.type(), pt2.type()));
    vlPairs.add(pt1, pt2); // add pair to Verlet List
  }

  /*-------------------------------------------------------------*/

  int VerletList::totalSize() const
  {
    System& system = getSystemRef();
    int size = localSize();
    int allsize;

    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
  }

  int VerletList::localSize() const
  {
    System& system = getSystemRef();
    // return vlPairs.size();
    return num_pairs;
  }

  python::tuple VerletList::getPair(int i) {
    if (i <= 0 || i > vlPairs.size()) {
      std::cout << "ERROR VerletList pair " << i << " does not exists" << std::endl;
      return python::make_tuple();
    } else {
      return python::make_tuple(vlPairs[i-1].first->id(), vlPairs[i-1].second->id());
    }
  }


  bool VerletList::exclude(longint pid1, longint pid2) {

      throw std::runtime_error("Exclusions in VerletList not implemented.");

      exList.insert(std::make_pair(pid1, pid2));

      return true;
  }


  /*-------------------------------------------------------------*/

  VerletList::~VerletList()
  {
    LOG4ESPP_INFO(theLogger, "~VerletList");

    if (!connectionResort.connected()) {
      connectionResort.disconnect();
    }
  }

  /*-------------------------------------------------------------*/

  void VerletList::resetTimers()
  {
    timeRebuild = 0.0;
  }

  void VerletList::loadTimers(real* t)
  {
    t[0] = timeRebuild;
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

  void VerletList::registerPython() {
    using namespace espressopp::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2)
          = &VerletList::exclude;


    class_<VerletList, shared_ptr<VerletList> >
      ("vectorization_VerletList", init< shared_ptr<System>, shared_ptr<Vectorization>, real, bool, int >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
      .def("totalSize", &VerletList::totalSize)
      .def("localSize", &VerletList::localSize)
      .def("getPair", &VerletList::getPair)
      .def("exclude", pyExclude)
      .def("rebuild", &VerletList::rebuild)
      .def("rebuildPairs", &VerletList::rebuildPairs)
      .def("connect", &VerletList::connect)
      .def("disconnect", &VerletList::disconnect)
      .def("getVerletCutoff", &VerletList::getVerletCutoff)
      .def("resetTimers", &VerletList::resetTimers)
      .def("getTimers", wrapGetTimers)
      ;
  }

}}
