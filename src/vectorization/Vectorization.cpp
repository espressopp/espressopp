/*
  Copyright (C) 2019
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

#include "Vectorization.hpp"
#include "storage/Storage.hpp"

namespace espressopp {
  namespace vectorization {

    using integrator::MDIntegrator;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    LOG4ESPP_LOGGER(Vectorization::logger, "Vectorization");

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// constructor
    Vectorization::Vectorization(
      shared_ptr<System> system,
      shared_ptr<MDIntegrator> mdintegrator,
      Mode mode
      ): SystemAccess(system), mdintegrator(mdintegrator), mode(mode)
    {
      std::string mode_str = (mode==ESPP_VEC_AOS) ? "AOS" : "SOA";
      LOG4ESPP_INFO(logger,"Using vectorization mode: " << mode_str);
      connect();
      resetCells(); // immediately retrieve cell information
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// destructor
    Vectorization::~Vectorization()
    {
      disconnect();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// connect to boost signals in integrator and storage
    void Vectorization::connect()
    {
      sigResetParticles     = getSystem()->storage->onParticlesChanged.connect(
                                boost::signals2::at_front, // call first due to reordering
                                boost::bind(&Vectorization::resetParticles, this));
      sigResetCells         = getSystem()->storage->onCellAdjust.connect(
                                boost::signals2::at_back,
                                boost::bind(&Vectorization::resetCells, this));
      sigBefCalcForces      = mdintegrator->aftInitF.connect(
                                boost::signals2::at_back,
                                boost::bind(&Vectorization::befCalcForces,this));
      sigUpdateForces       = mdintegrator->aftCalcFLocal.connect(
                                boost::signals2::at_front,
                                boost::bind(&Vectorization::updateForces,this));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// disconnect boost signals made with connect()
    void Vectorization::disconnect()
    {
      sigResetParticles.disconnect();
      sigResetCells.disconnect();
      sigBefCalcForces.disconnect();
      sigUpdateForces.disconnect();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// reset and update particleArray from current storage
    void Vectorization::resetParticles()
    {
      particleArray.copy_from(getSystem()->storage->getLocalCells(), mode);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// reset and update cell mapping and neighbor lists from current storage
    void Vectorization::resetCells()
    {
      neighborList = CellNeighborList(getSystem()->storage);
      LOG4ESPP_TRACE(logger,"neighborList, ncells: "<<neighborList.numCells()<<" nnbrs: "<<neighborList.maxNumNeighbors());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// set force array/s to zero
    void Vectorization::befCalcForces()
    {
      if(mode==ESPP_VEC_AOS)
      {
        auto& f = particleArray.force;
        real* el = &(f[0].x);
        size_t end = 4*f.size();

        #pragma vector always
        #pragma vector aligned
        #pragma ivdep
        for(size_t i=0; i<end; i++) el[i] = 0.0;
      }
      else
      {
        auto& f_x = particleArray.f_x;
        auto& f_y = particleArray.f_y;
        auto& f_z = particleArray.f_z;
        std::fill(f_x.begin(),f_x.end(),0.0);
        std::fill(f_y.begin(),f_y.end(),0.0);
        std::fill(f_z.begin(),f_z.end(),0.0);
      }

      // overwrite particleArray positon data
      particleArray.update_from_position_only(getSystem()->storage->getLocalCells());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// add forces back to storage
    void Vectorization::updateForces()
    {
      // add particle array forces back to localCells
      particleArray.add_to_force_only(getSystem()->storage->getLocalCells());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// registration with python
    void Vectorization::registerPython()
    {
      using namespace espressopp::python;

      class_<Vectorization, shared_ptr<Vectorization> >
        ("Vectorization", init< shared_ptr<System>, shared_ptr<MDIntegrator>, Mode >())
        .def(init< shared_ptr<System>, shared_ptr<MDIntegrator> >())
        ;

      enum_<Mode>("VectorizationMode")
        .value("SOA",ESPP_VEC_SOA)
        .value("AOS",ESPP_VEC_AOS)
        ;
    }

  }
}
