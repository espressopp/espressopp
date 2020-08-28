/*
  Copyright (C) 2019-2020
      Max Planck Institute for Polymer Research & JGU Mainz

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

#include "ParticleArray.hpp"
#include "Cell.hpp"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////////////////////////
namespace espressopp { namespace vectorization
{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::copyFrom(CellList const& srcCells, Mode mode_)
  {
    mode = mode_;

    size_t const numCells = srcCells.size();
    cellRange_.clear();
    cellRange_.reserve(numCells+1);
    sizes_.clear();
    sizes_.reserve(numCells);

    // preallocate
    size_t total_size = 0, total_data_size = 0;

    for(size_t ic=0; ic<numCells; ic++)
    {
      Cell* cell = srcCells[ic];
      cellRange_.push_back(total_data_size);
      size_t const size = cell->particles.size();
      sizes_.push_back(size);
      // require only at least one particle for padding
      size_t const data_size = (size+1);

      total_size += size;
      total_data_size += data_size;
    }
    cellRange_.push_back(total_data_size);

    if(total_data_size>reserve_size_)
    {
      reserve_size_ = 2*total_data_size;
      if(mode==ESPP_VEC_AOS)
      {
        ESPP_PARTICLEARRAY_AOS_APPLY(resize(reserve_size_));
      }
      else if(mode==ESPP_VEC_SOA)
      {
        ESPP_PARTICLEARRAY_SOA_APPLY(resize(reserve_size_));
      }
      else
      {
        throw std::runtime_error("ParticleArray: Given mode not implemented.");
      }
      // std::cout << __FUNCTION__ << ": Reallocated " << total_data_size << std::endl;
      // TODO: implement as LOG4ESPP_WARN
    }

    size_ = total_size;
    data_size_ = total_data_size;

    // fill values
    for(size_t ic=0; ic<numCells; ic++)
    {
      size_t start = cellRange_[ic];
      ParticleList const& particlelist = srcCells[ic]->particles;
      updateFrom(particlelist, start);
    }

    if(mode==ESPP_VEC_AOS)
    {
      for(size_t ic=0; ic<numCells; ic++)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        size_t end = cellRange_[ic] + particlelist.size();
        size_t data_end = cellRange_[ic+1];
        for(size_t ip=end; ip<data_end; ip++) position[ip]={large_pos,large_pos,large_pos,0};
      }
    }
    else if(mode==ESPP_VEC_SOA)
    {
      for(size_t ic=0; ic<numCells; ic++)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        size_t end = cellRange_[ic] + particlelist.size();
        size_t data_end = cellRange_[ic+1];
        for(size_t ip=end; ip<data_end; ip++) p_x  [ip] = large_pos;
        for(size_t ip=end; ip<data_end; ip++) p_y  [ip] = large_pos;
        for(size_t ip=end; ip<data_end; ip++) p_z  [ip] = large_pos;
        for(size_t ip=end; ip<data_end; ip++) type [ip] = 0;
      }
    }
    else
    {
      throw std::runtime_error("ParticleArray: Given mode not implemented.");
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::updateFrom(std::vector<Particle> const& particlelist, size_t start)
  {
    size_t end = start + particlelist.size();

    if(mode==ESPP_VEC_AOS)
    {
      for(size_t pi=start, pli=0; pi<end; pi++,pli++)
      {
        Particle const& p = particlelist[pli];
        position [pi] = Real3DInt(p.position(),p.type());
        force    [pi] = p.force();
      }
    }
    else
    {
      for(size_t pi=start, pli=0; pi<end; pi++,pli++)
      {
        Particle const& p = particlelist[pli];
        type [pi] = p.type();
        p_x  [pi] = p.position()[0];
        p_y  [pi] = p.position()[1];
        p_z  [pi] = p.position()[2];
        f_x  [pi] = p.force()[0];
        f_y  [pi] = p.force()[1];
        f_z  [pi] = p.force()[2];
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::updateFromPositionOnly(CellList const& srcCells)
  {
    if(mode==ESPP_VEC_AOS)
    {
      for(size_t ic=0; ic<srcCells.size(); ic++)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        const size_t start = cellRange_[ic];
        const size_t end = start + particlelist.size();
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          position[pi] = particlelist[pli].position();
        }
      }
    }
    else
    {
      for(size_t ic=0; ic<srcCells.size(); ic++)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        const size_t start = cellRange_[ic];
        const size_t end = start + particlelist.size();
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          Particle const& p = particlelist[pli];
          p_x[pi] = p.position()[0];
          p_y[pi] = p.position()[1];
          p_z[pi] = p.position()[2];
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::addToForceOnly(CellList & srcCells) const
  {
    if(mode==ESPP_VEC_AOS)
    {
      for(size_t ic=0; ic<srcCells.size(); ic++)
      {
        ParticleList & particlelist = srcCells[ic]->particles;
        const size_t start = cellRange_[ic];
        const size_t end = start + particlelist.size();
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].force() += force[pi].to_Real3D();
        }
      }
    }
    else
    {
      for(size_t ic=0; ic<srcCells.size(); ic++)
      {
        ParticleList & particlelist = srcCells[ic]->particles;
        const size_t start = cellRange_[ic];
        const size_t end = start + particlelist.size();
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].force() += Real3D(f_x[pi],f_y[pi],f_z[pi]);
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  bool ParticleArray::checkSizes()
  {
    if(size_>data_size_) return false;

    if(mode==ESPP_VEC_AOS)
    {
      if(position    .size()<data_size_) return false;
      if(force       .size()<data_size_) return false;
    }
    else
    {
      if(p_x         .size()<data_size_) return false;
      if(p_y         .size()<data_size_) return false;
      if(p_z         .size()<data_size_) return false;
      if(f_x         .size()<data_size_) return false;
      if(f_y         .size()<data_size_) return false;
      if(f_z         .size()<data_size_) return false;
      if(type        .size()<data_size_) return false;
    }

    return true;
  }

}}
