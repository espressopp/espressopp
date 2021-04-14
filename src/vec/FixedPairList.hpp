/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013,2016
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

// ESPP_CLASS
#ifndef VEC_FIXEDPAIRLIST_HPP
#define VEC_FIXEDPAIRLIST_HPP

#include "Vectorization.hpp"

#include "log4espp.hpp"
#include "python.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espressopp {
  namespace vec {

    typedef std::pair<size_t,size_t> ParticlePair;
    typedef AlignedVector<ParticlePair> PairList;

    class FixedPairList
      : public PairList
    {
    public:
      typedef boost::unordered_multimap<size_t, size_t> GlobalPairs;
      shared_ptr <Vectorization> vectorization;

    protected:
      boost::signals2::connection sigBeforeSend;
      boost::signals2::connection sigOnLoadCells;
      boost::signals2::connection sigAfterRecv;
      GlobalPairs globalPairs;
      real longtimeMaxBondSqr;

    public:
      FixedPairList(shared_ptr<Vectorization> vectorization);
      ~FixedPairList();

      real getLongtimeMaxBondSqr();
      void setLongtimeMaxBondSqr(real d);
      void resetLongtimeMaxBondSqr();

      /** Add the given particle pair to the list on this processor if the
          particle with the lower id belongs to this processor.  Note that
          this routine does not check whether the pair is inserted on
          another processor as well.
          \return whether the particle was inserted on this processor.
      */
      bool add(size_t pid1, size_t pid2);
      void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
      void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
      void onLoadCells();
      void remove();
      std::vector<size_t> getPairList();
      python::list getBonds();
      python::list getAllBonds();
      GlobalPairs* getGlobalPairs() {return &globalPairs;};

      /** Get the number of bonds in the GlobalPairs list */
      int size() {
          return globalPairs.size();
      }

      int totalSize();

      static void registerPython();

    private:
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };

  }
}

#endif//VEC_FIXEDPAIRLIST_HPP
