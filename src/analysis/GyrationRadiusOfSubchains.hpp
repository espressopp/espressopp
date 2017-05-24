/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research
  
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
#ifndef _ANALYSIS_GYRATIONRADIUS_HPP
#define _ANALYSIS_GYRATIONRADIUS_HPP

#include "ConfigsParticleDecomp.hpp"

namespace espressopp {
  namespace analysis {

    /*
     * Class derives from ConfigsParticleDecomp
    */
    class GyrationRadiusOfSubchains : public ConfigsParticleDecomp {

    public:
      
      GyrationRadiusOfSubchains(shared_ptr<System> system, int _chainlength): ConfigsParticleDecomp (system, _chainlength){
        // by default calculation progress is printed
        //setPrint_progress(true);
        setPrint_progress(false);

	//key = "unfolded";

        int localN = system -> storage -> getNRealParticles();
	//cout << "The number of real particles " << localN << endl;
        boost::mpi::all_reduce(*system->comm, localN, num_of_part, std::plus<int>());
	//cout << "Total of real particles " << num_of_part << endl;

	int localNN = system -> storage -> getNLocalParticles();
	int num_of_local;
	//cout << "The number of local particles " << localNN << endl;
	boost::mpi::all_reduce(*system->comm, localNN, num_of_local, std::plus<int>());
	//cout << "Total of local particles " << num_of_local << endl;
        
        int n_nodes = system -> comm -> size();
        int this_node = system -> comm -> rank();
              
        //for monodisperse chains
        int num_chains = num_of_part / chainlength;
        int local_num_chains = (int) ceil( (double)num_chains / n_nodes );

        int local_num_of_part = num_of_part / n_nodes + 1;

	idToCpu.clear();
        int nodeNum = 0;
        int count = 0;
	for(long unsigned int k = 1; k <= num_chains ;k++){
          idToCpu[k] = nodeNum;
          count ++;
          if(count>=local_num_of_part){
            count = 0;
            nodeNum++;
          }
        }
      }
      ~GyrationRadiusOfSubchains() {}
      
      virtual python::list compute() const;

      void setPrint_progress(bool _print_progress){
        print_progress = _print_progress;
      }
      bool getPrint_progress(){return print_progress;}
      
      static void registerPython();
    private:
      bool print_progress;
    };
  }
}

#endif
