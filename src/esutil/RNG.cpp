/*
  Copyright (C) 2012,2013,2019
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
#include "RNG.hpp"
#include "mpi.hpp"
#include "types.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "boost/serialization/string.hpp"

using namespace boost;

namespace espressopp {
  namespace esutil {

    RNG::RNG(long _seed): boostRNG(make_shared< RNGType >(_seed + mpiWorld->rank())),
            normalVariate(*boostRNG, normal_distribution< real >(0.0, 1.0)),
            uniformOnSphereVariate(*boostRNG, uniform_on_sphere< real, Real3D >(3)),
            seed_(_seed)
    		//gammaVariate(*boostRNG, gamma_distribution< real >(1, 1.0)), //TODO this line is nonsense: alpha=1 is trivial
    {}

    void RNG::seed(long _seed) {
      // Seed the RNG for the given CPU
      boostRNG->seed(_seed + mpiWorld->rank());
      seed_ = _seed;
    }

    long RNG::get_seed() {
      return seed_;
    }

    real RNG::operator()() {
      variate_generator< RNGType&, uniform_01<> > uni(*boostRNG, uniform_01<>());
      return uni();
    }

    int RNG::operator()(int N) {
      uniform_smallint< int > uni_dist(0, N-1);
      variate_generator< RNGType&, uniform_smallint< int > > uni(*boostRNG, uni_dist);
      return uni();
    }

    real RNG::normal() {
      return normalVariate();
    }

    real RNG::gamma(unsigned int alpha) {
    	gamma_distribution< real > gamma_dist(alpha, 1.0); //scale parameter \beta=1.0
    	variate_generator< RNGType&, gamma_distribution< real > > gamma_var(*boostRNG, gamma_dist);
    	return gamma_var();
    }


    Real3D RNG::uniformOnSphere() {
      return uniformOnSphereVariate();
    }

    shared_ptr< RNGType > RNG::getBoostRNG() {
      return boostRNG;
    }

    void RNG::saveState(long long step) {
      std::ostringstream oss;
      oss << *boostRNG;
      if ( mpiWorld->rank()==0 ) {
        std::vector<std::string> oss_str_vec;
        boost::mpi::gather(*mpiWorld, oss.str(), oss_str_vec, 0);
        std::string outfn = "rng." + std::to_string(step);
        std::ofstream ofs(outfn, std::ios::out | std::ios::binary);
        if (ofs.is_open()) {
          for(unsigned int i = 0;i<oss_str_vec.size();++i) {
            ofs << oss_str_vec[i] << std::endl;
          }
        }
        ofs.close();
      } else {
        boost::mpi::gather(*mpiWorld, oss.str(), 0);
      }
    }

    void RNG::loadState(const char* fname) {
      std::ifstream ifs;
      std::string rng_str;
      std::vector<std::string> rng_str_vec;
      unsigned int rng_size;
      if ( mpiWorld->rank()==0 ) {
        ifs.open(fname, std::ios::in | std::ios::binary);
        std::string token;
        while(std::getline(ifs, token)) {
          rng_str_vec.push_back(token);
          token.clear();
        }
        rng_size = rng_str_vec.size();
        if (ifs.is_open()) ifs.close();
      }
      boost::mpi::broadcast(*mpiWorld, rng_size, 0);
      if (rng_size != static_cast<unsigned int>(mpiWorld->size())) {
        std::cerr << "# Warning: The number of processes does not match the number of saved random number generators."  << std::endl << "           The old states are not restored.\n";
        return;
      }
      boost::mpi::scatter(*mpiWorld, rng_str_vec, rng_str, 0);
      std::istringstream iss(rng_str);
      iss >> *boostRNG;
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    RNG::registerPython() {
      using namespace espressopp::python;

      real (RNG::*pyCall1)() = &RNG::operator();
      int (RNG::*pyCall2)(int) = &RNG::operator();


      class_< RNG >("esutil_RNG", init< boost::python::optional< long > >())
        .def("seed", &RNG::seed)
        .def("__call__", pyCall1)
        .def("__call__", pyCall2)
        .def("normal", &RNG::normal)
        .def("gamma", &RNG::gammaOf1)
        .def("gamma", &RNG::gamma)
        .def("uniformOnSphere", &RNG::uniformOnSphere)
        .def("get_seed", &RNG::get_seed)
        .def("saveState", &RNG::saveState)
        .def("loadState", &RNG::loadState);
    }
  }
}
