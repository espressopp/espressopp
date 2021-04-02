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

#ifndef VEC_STORAGEVEC_HPP
#define VEC_STORAGEVEC_HPP

#include "vec/include/types.hpp"

#include "python.hpp"
#include "log4espp.hpp"

#include <boost/signals2.hpp>

namespace espressopp { namespace vec {
  namespace storage {

    class StorageVec
    {
    public:
      StorageVec(shared_ptr<Vectorization> vectorization);

      virtual void loadCells() = 0;

      virtual void unloadCells() = 0;

      boost::signals2::signal<void ()> onLoadCells;

      static void registerPython();

    protected:
      shared_ptr<Vectorization> vectorization;

    private:
      static LOG4ESPP_DECL_LOGGER(logger);
    };

  }
}}

#endif//VEC_STORAGEVEC_HPP
