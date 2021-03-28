/*
  Copyright (C) 2021
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

#define ESPP_PARTICLEARRAY_AOS_APPLY(COMMAND) \
  position    . COMMAND ; \
  velocity    . COMMAND ; \
  force       . COMMAND ; \
  id          . COMMAND ; \
  mass        . COMMAND ; \
  q           . COMMAND ; \
  /* */

#define ESPP_PARTICLEARRAY_SOA_APPLY(COMMAND) \
  p_x         . COMMAND ; \
  p_y         . COMMAND ; \
  p_z         . COMMAND ; \
  v_x         . COMMAND ; \
  v_y         . COMMAND ; \
  v_z         . COMMAND ; \
  f_x         . COMMAND ; \
  f_y         . COMMAND ; \
  f_z         . COMMAND ; \
  type        . COMMAND ; \
  id          . COMMAND ; \
  mass        . COMMAND ; \
  q           . COMMAND ; \
  /* */

///////////////////////////////////////////////////////////////////////////////////////////////////
namespace espressopp { namespace vec
{
  ParticleArray::ParticleArray()
  {}

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::markRealCells(CellList const& realcellsIn, const Cell* cell0)
  {
    realCells_.clear();
    realCells_.reserve(realcellsIn.size());
    for(const Cell* cell: realcellsIn)
    {
      std::ptrdiff_t icell = cell - cell0;
      if(icell<0) {
        std::cerr << "ParticleArray::markRealCells " <<
          "Provided cell0 must be less than address of realCells" << std::endl;
      }
      realCells_.push_back(icell);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::markRealCells(std::vector<size_t> const& realcellsIn)
  {
    realCells_ = realcellsIn;
  }

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
      // // require only at one particle for padding
      // size_t const data_size = (size+1);
      // require only at least one particle for padding but fill one full row
      size_t const data_size = calc_data_size(size+1);

      total_size += size;
      total_data_size += data_size;
    }
    cellRange_.push_back(total_data_size);

    if(total_data_size>reserve_size_)
    {
      reserve_size_ = 2*total_data_size;
      if(mode==ESPP_VEC_AOS)
      {
        ESPP_PARTICLEARRAY_AOS_APPLY(clear());
        ESPP_PARTICLEARRAY_AOS_APPLY(resize(reserve_size_));
      }
      else if(mode==ESPP_VEC_SOA)
      {
        ESPP_PARTICLEARRAY_SOA_APPLY(clear());
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
    auto updateCell = [this,&srcCells](size_t ic)
    {
      size_t start = cellRange_[ic];
      ParticleList const& particlelist = srcCells[ic]->particles;
      updateFrom(particlelist, start);
    };
    for(size_t ic=0; ic<numCells; ic++) updateCell(ic);

    if(mode==ESPP_VEC_AOS)
    {
      auto fillCell = [this,&srcCells](size_t ic)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        const size_t end = cellRange_[ic] + particlelist.size();
        const size_t data_end = cellRange_[ic+1];
        for(size_t ip=end; ip<data_end; ip++) position[ip]={large_pos,large_pos,large_pos,0};
      };
      for(size_t ic=0; ic<numCells; ic++) fillCell(ic);
    }
    else if(mode==ESPP_VEC_SOA)
    {
      auto fillCell = [this,&srcCells](size_t ic)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        const size_t end = cellRange_[ic] + particlelist.size();
        const size_t data_end = cellRange_[ic+1];
        for(size_t ip=end; ip<data_end; ip++) p_x  [ip] = large_pos;
        for(size_t ip=end; ip<data_end; ip++) p_y  [ip] = large_pos;
        for(size_t ip=end; ip<data_end; ip++) p_z  [ip] = large_pos;
        for(size_t ip=end; ip<data_end; ip++) type [ip] = 0;
      };
      for(size_t ic=0; ic<numCells; ic++) fillCell(ic);
    }
    else
    {
      throw std::runtime_error("ParticleArray: Given mode not implemented.");
    }

    {
      auto fillCell = [this,&srcCells](size_t ic)
      {
        ParticleList const& particlelist = srcCells[ic]->particles;
        size_t end = cellRange_[ic] + particlelist.size();
        size_t data_end = cellRange_[ic+1];
        for(size_t ip=end; ip<data_end; ip++) id    [ip] = -1;
        for(size_t ip=end; ip<data_end; ip++) mass  [ip] = 1.0;
        for(size_t ip=end; ip<data_end; ip++) q     [ip] = 0.0;
      };
      for(size_t ic=0; ic<numCells; ic++) fillCell(ic);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::updateFromPositionVelocity(CellList const& srcCells, bool realOnly)
  {
    const size_t numCells = srcCells.size();

    auto updateCell = [this,&srcCells](size_t ic)
    {
      ParticleList const& particlelist = srcCells[ic]->particles;
      const size_t start = cellRange_[ic];
      const size_t end = start + particlelist.size();
      if(mode==ESPP_VEC_AOS)
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          position[pi] = particlelist[pli].position();
          velocity[pi] = particlelist[pli].velocity();
        }
      }
      else
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          Particle const& p = particlelist[pli];
          p_x[pi] = p.position()[0];
          p_y[pi] = p.position()[1];
          p_z[pi] = p.position()[2];
          v_x[pi] = p.velocity()[0];
          v_y[pi] = p.velocity()[1];
          v_z[pi] = p.velocity()[2];
        }
      }
    };

    if(realOnly)
    {
      auto updateRealCell = [this,&updateCell](size_t ir)
      {
        updateCell(realCells_.at(ir));
      };
      const size_t numRealCells = realCells_.size();
      for(size_t ir=0; ir<numRealCells; ir++) updateRealCell(ir);
    }
    else
    {
      for(size_t ic=0; ic<numCells; ic++) updateCell(ic);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::updateToPositionVelocity(CellList & srcCells, bool realOnly) const
  {
    /// copy particle info (position and velocity ) back to srcCells
    /// TODO: track which properties were modified using flags and possibly implement templates to
    /// offload each modified property

    /// TODO: enable option to update only real cells since ghost cells will be cleared during rebuild

    /// NOTE: Currently the only properties that needs to update back is the position and velocity;
    /// forces are modified but are not used after the update, while type, id, mass and q are
    /// considered static throughout the simulation.

    /// Verify source cell sizes
    /// Deactivate in production builds
    #if defined(ESPP_VEC_DEBUG)
      verify(srcCells);
    #endif

    /// TODO: Remove block later
    #if 0
    if(mode==ESPP_VEC_AOS)
    {
      auto updateCell = [this,&srcCells](size_t ic)
      {
        ParticleList & particlelist = srcCells[ic]->particles;
        const size_t start = cellRange_[ic];
        const size_t end = start + particlelist.size();
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].position() = position[pi].to_Real3D();
          particlelist[pli].velocity() = velocity[pi].to_Real3D();
        }
      };
      for(size_t ic=0; ic<numCells; ic++) updateCell(ic);
    }
    else
    {
      auto updateCell = [this,&srcCells](size_t ic)
      {
        ParticleList & particlelist = srcCells[ic]->particles;
        const size_t start = cellRange_[ic];
        const size_t end = start + particlelist.size();
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].position() = Real3D(p_x[pi],p_y[pi],p_z[pi]);
          particlelist[pli].velocity() = Real3D(v_x[pi],v_y[pi],v_z[pi]);
        }
      };
      for(size_t ic=0; ic<numCells; ic++) updateCell(ic);
    }
    #endif

    auto updateCell = [this,&srcCells](size_t ic)
    {
      ParticleList & particlelist = srcCells[ic]->particles;
      const size_t start = cellRange_[ic];
      const size_t end = start + particlelist.size();
      if(mode==ESPP_VEC_AOS)
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].position() = position[pi].to_Real3D();
          particlelist[pli].velocity() = velocity[pi].to_Real3D();
        }
      }
      else
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].position() = Real3D(p_x[pi],p_y[pi],p_z[pi]);
          particlelist[pli].velocity() = Real3D(v_x[pi],v_y[pi],v_z[pi]);
        }
      }
    };

    if(realOnly)
    {
      auto updateRealCell = [this,&updateCell](size_t ir)
      {
        updateCell(realCells_.at(ir));
      };
      const size_t numRealCells = realCells_.size();
      for(size_t ir=0; ir<numRealCells; ir++) updateRealCell(ir);
    }
    else
    {
      const size_t numCells = srcCells.size();
      for(size_t ic=0; ic<numCells; ic++) updateCell(ic);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::verify(CellList const& srcCells) const
  {
    if((cellRange_.size() != srcCells.size()+1) || (sizes_.size() != srcCells.size())) {
      std::ostringstream oss;
      oss << "ParticleArray::verify " <<
        "Incompatible number of cells in srcCells";
      throw std::runtime_error(oss.str());
    }
    if(!checkSizes()) {
      std::ostringstream oss;
      oss << "ParticleArray::verify " <<
        "Incompatible number of entries in ParticleArray found with checkSizes";
      throw std::runtime_error(oss.str());
    }
    for(size_t ic=0; ic<srcCells.size(); ic++)
    {
      ParticleList & particlelist = srcCells[ic]->particles;
      if(sizes_[ic] != particlelist.size()) {
        std::ostringstream oss;
        oss << "ParticleArray::verify " <<
          "Incompatible number of particles in srcCells[" << ic << "]: " <<
          "expected "<< sizes_[ic] <<", got "<<particlelist.size();
        throw std::runtime_error(oss.str());
      }
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

        id       [pi] = p.id();
        mass     [pi] = p.mass();
        q        [pi] = p.q();

        position [pi] = Real3DInt(p.position(),p.type());
        velocity [pi] = p.velocity();
        force    [pi] = p.force();
      }
    }
    else
    {
      for(size_t pi=start, pli=0; pi<end; pi++,pli++)
      {
        Particle const& p = particlelist[pli];

        id    [pi] = p.id();
        mass  [pi] = p.mass();
        q     [pi] = p.q();

        type  [pi] = p.type();
        p_x   [pi] = p.position()[0];
        p_y   [pi] = p.position()[1];
        p_z   [pi] = p.position()[2];

        v_x   [pi] = p.velocity()[0];
        v_y   [pi] = p.velocity()[1];
        v_z   [pi] = p.velocity()[2];

        f_x   [pi] = p.force()[0];
        f_y   [pi] = p.force()[1];
        f_z   [pi] = p.force()[2];
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
  void ParticleArray::updateToPositionOnly(CellList & srcCells) const
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
          particlelist[pli].position() = position[pi].to_Real3D();
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
          particlelist[pli].position() = Real3D(p_x[pi],p_y[pi],p_z[pi]);
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::updateFromForceOnly(CellList const& srcCells, bool realOnly)
  {
    auto updateCell = [this,&srcCells](size_t ic)
    {
      ParticleList const& particlelist = srcCells[ic]->particles;
      const size_t start = cellRange_[ic];
      const size_t end = start + particlelist.size();
      if(mode==ESPP_VEC_AOS)
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          force[pi] = particlelist[pli].force();
        }
      }
      else
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          Particle const& p = particlelist[pli];
          f_x[pi] = p.force()[0];
          f_y[pi] = p.force()[1];
          f_z[pi] = p.force()[2];
        }
      }
    };

    if(realOnly)
    {
      auto updateRealCell = [this,&updateCell](size_t ir)
      {
        updateCell(realCells_.at(ir));
      };
      const size_t numRealCells = realCells_.size();
      for(size_t ir=0; ir<numRealCells; ir++) updateRealCell(ir);
    }
    else
    {
      const size_t numCells = srcCells.size();
      for(size_t ic=0; ic<numCells; ic++) updateCell(ic);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticleArray::updateToForceOnly(CellList & srcCells, bool realOnly) const
  {
    auto updateCell = [this,&srcCells](size_t ic)
    {
      ParticleList & particlelist = srcCells[ic]->particles;
      const size_t start = cellRange_[ic];
      const size_t end = start + particlelist.size();
      if(mode==ESPP_VEC_AOS)
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].force() = force[pi].to_Real3D();
        }
      }
      else
      {
        for(size_t pi=start, pli=0; pi<end; pi++,pli++)
        {
          particlelist[pli].force() = Real3D(f_x[pi],f_y[pi],f_z[pi]);
        }
      }
    };

    if(realOnly)
    {
      auto updateRealCell = [this,&updateCell](size_t ir)
      {
        updateCell(realCells_.at(ir));
      };
      const size_t numRealCells = realCells_.size();
      for(size_t ir=0; ir<numRealCells; ir++) updateRealCell(ir);
    }
    else
    {
      const size_t numCells = srcCells.size();
      for(size_t ic=0; ic<numCells; ic++) updateCell(ic);
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
  bool ParticleArray::checkSizes() const
  {
    if(size_>data_size_) return false;

    if(id   .size()<data_size_) return false;
    if(mass .size()<data_size_) return false;
    if(q    .size()<data_size_) return false;

    if(mode==ESPP_VEC_AOS)
    {
      if(position .size()<data_size_) return false;
      if(velocity .size()<data_size_) return false;
      if(force    .size()<data_size_) return false;
    }
    else
    {
      if(p_x      .size()<data_size_) return false;
      if(p_y      .size()<data_size_) return false;
      if(p_z      .size()<data_size_) return false;
      if(v_x      .size()<data_size_) return false;
      if(v_y      .size()<data_size_) return false;
      if(v_z      .size()<data_size_) return false;
      if(f_x      .size()<data_size_) return false;
      if(f_y      .size()<data_size_) return false;
      if(f_z      .size()<data_size_) return false;
      if(type     .size()<data_size_) return false;
    }

    return true;
  }

}}
