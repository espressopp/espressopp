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

#include "Vectorization.hpp"
#include "storage/Storage.hpp"
#include "vec/integrator/MDIntegratorVec.hpp"

namespace espressopp {
  namespace vec {

    ///////////////////////////////////////////////////////////////////////////////////////////////
    LOG4ESPP_LOGGER(Vectorization::logger, "Vectorization");

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// constructor
    Vectorization::Vectorization(
      shared_ptr<System> _system,
      shared_ptr<MDIntegrator> _mdintegrator,
      Mode _vecMode
      )
      : SystemAccess(_system)
      , mdintegrator(_mdintegrator)
      , vecMode(_vecMode)
      , vecLevel(1)
    {
      std::string mode_str = (vecMode==ESPP_VEC_AOS) ? "AOS" : "SOA";
      LOG4ESPP_INFO(logger,"Using vectorization mode: " << mode_str << " level " << vecLevel);
      connect_level_1();
      resetCells(); // immediately retrieve cell information
      std::cout << "Vectorization level " << vecLevel << std::endl;
    }

    Vectorization::Vectorization(
      shared_ptr<System> _system,
      Mode _vecMode
      )
      : SystemAccess(_system)
      , vecMode(_vecMode)
      , vecLevel(2)
    {
      std::string mode_str = (vecMode==ESPP_VEC_AOS) ? "AOS" : "SOA";
      LOG4ESPP_INFO(logger,"Using vectorization mode: " << mode_str << " level " << vecLevel);
      connect_level_2();
      // resetCells(); // TODO: retrieve cell information when connecting to Storage
      std::cout << "Vectorization level " << vecLevel << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// destructor
    Vectorization::~Vectorization()
    {
      disconnect();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// connect to boost signals in integrator and storage
    void Vectorization::connect_level_1()
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

    void Vectorization::connect_level_2()
    {
      // throw std::runtime_error("Vectorization level>=2 not implemented.");
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
    /// reset and update particles from current storage
    void Vectorization::resetParticles()
    {
      particles.copyFrom(getSystem()->storage->getLocalCells(), vecMode);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// reset and update cell mapping and neighbor lists from current storage
    void Vectorization::resetCells()
    {
      if(!getSystem()->storage){
        throw std::runtime_error("System has no storage");
      }
      resetCells(getSystem()->storage.get());
    }

    void Vectorization::resetCells(Storage* storage)
    {
      CellList const& localCells = storage->getLocalCells();
      CellList const& realCells  = storage->getRealCells();
      Cell* const cell0 = localCells[0];
      std::vector<size_t> realCellIdx;
      for(Cell* rc: realCells) {
        realCellIdx.push_back(rc-cell0);
      }
      particles.markRealCells(realCellIdx);
      neighborList = CellNeighborList(cell0, localCells, realCellIdx);
      LOG4ESPP_TRACE(logger,"neighborList, ncells: "<<neighborList.numCells()
        <<" nnbrs: "<<neighborList.maxNumNeighbors());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// set force array/s to zero
    void Vectorization::befCalcForces()
    {
      if(vecMode==ESPP_VEC_AOS)
      {
        auto& f = particles.force;
        real* el = &(f[0].x);
        size_t end = 4*f.size();

        #pragma vector always
        #pragma vector aligned
        #pragma ivdep
        for(size_t i=0; i<end; i++) el[i] = 0.0;
      }
      else
      {
        auto& f_x = particles.f_x;
        auto& f_y = particles.f_y;
        auto& f_z = particles.f_z;
        std::fill(f_x.begin(),f_x.end(),0.0);
        std::fill(f_y.begin(),f_y.end(),0.0);
        std::fill(f_z.begin(),f_z.end(),0.0);
      }

      // overwrite particles position data
      particles.updateFromPositionOnly(getSystem()->storage->getLocalCells());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// add forces back to storage
    void Vectorization::updateForces()
    {
      // add particle array forces back to localCells
      particles.addToForceOnly(getSystem()->storage->getLocalCells());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// registration with python
    void Vectorization::registerPython()
    {
      using namespace espressopp::python;

      class_<Vectorization, shared_ptr<Vectorization> >
        ("vec_Vectorization",
             init< shared_ptr<System>, shared_ptr<MDIntegrator>, Mode >())
        .def(init< shared_ptr<System>, shared_ptr<MDIntegrator> >())
        .def(init< shared_ptr<System>, Mode >())
        .def(init< shared_ptr<System>>())
        .add_property("level", &Vectorization::getVecLevel)
        .def_readwrite("storageVec", &Vectorization::storageVec)
        ;

      enum_<Mode>("VecMode")
        .value("SOA",ESPP_VEC_SOA)
        .value("AOS",ESPP_VEC_AOS)
        ;
    }

  }
}
