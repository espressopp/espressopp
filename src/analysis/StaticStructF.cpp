/*
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
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "Configuration.hpp"
#include "StaticStructF.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"

#include <boost/serialization/map.hpp>

#include <math.h>       // cos and ceil and sqrt
#include <algorithm>    // std::min
#include <functional>   // std::plus
#include <time.h>       // time_t, for particle-distribution-to-cpu time

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

using namespace espressopp;
using namespace espressopp::iterator;
using namespace std;

namespace espressopp {
    namespace analysis {
        // currently only works for particles numbered like 0, 1, 2,...

        // nqx is a number which corresponds to the different x-values of the
        // diffraction vector q. greater nqx produces more different x-values  
        // bin_factor determines the size for the binning of q-vectors in using 
        // dq = 2*PI/boxlength as a reference value such that 
        // bin_size = bin_factor * dq 
        // more in detail:
        // dq is the shortest step of dqx, dqy, dqz - corresponding to the
        // longest side of the box. dq = min(dqx, dqy, dqz)
        // dqx, dqy, dqz are the cell length of the grid of possible q-vectors
        // dqx = 2*PI/Lx, dqy = 2*PI/Ly, dqz = 2*PI/Lz

        python::list StaticStructF::computeArray(int nqx, int nqy, int nqz,
                real bin_factor) const {
            time_t start;
            time(&start);
            cout << "collective calc starts " << ctime(&start) << "\n";
            //fist the system coords are saved at each CPU
            System& system = getSystemRef();
            esutil::Error err(system.comm);
            Real3D Li = system.bc->getBoxL(); //Box size (Lx, Ly, Lz)

            int nprocs = system.comm->size(); // number of CPUs
            int myrank = system.comm->rank(); // current CPU's number

            if (myrank == 0) {
                cout << "collective calc starts " << ctime(&start) << "\n";
            }

            int num_part = 0;
            ConfigurationPtr config = make_shared<Configuration > ();
            // loop over all CPU-numbers - to give all CPUs all particle coords
            for (int rank_i = 0; rank_i < nprocs; rank_i++) {
                map< size_t, Real3D > conf;
                if (rank_i == myrank) {
                    CellList realCells = system.storage->getRealCells();
                    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
                        int id = cit->id();
                        conf[id] = cit->position();
                    }
                }
                boost::mpi::broadcast(*system.comm, conf, rank_i);

                // for simplicity we will number the particles from 0
                for (map<size_t, Real3D>::iterator itr = conf.begin(); itr != conf.end(); ++itr) {
                    size_t id = itr->first;
                    Real3D p = itr->second;
                    config->set(id, p[0], p[1], p[2]);
                    //config->set(num_part, p[0], p[1], p[2]);
                    num_part++;
                }
            }
            if (myrank == 0) {
                time_t distributed;
                time(&distributed);
                cout << "particles on all CPUs " << ctime(&distributed) << "\n";
                cout << "distribution to CPUs took "
                        << difftime(distributed, start) << " seconds \n";
            }
            // now all CPUs have all particle coords and num_part is the total number
            // of particles

            // use all CPUs
            // TODO it could be a problem if   n_nodes > num_part

            // here starts calculation of the static structure factor

            //step size for qx, qy, qz
            real dqs[3];
            dqs[0] = 2. * M_PIl / Li[0];
            dqs[1] = 2. * M_PIl / Li[1];
            dqs[2] = 2. * M_PIl / Li[2];

            Real3D q;

            //calculations for binning
            real bin_size = bin_factor * min(dqs[0], (dqs[1], dqs[2]));
            real q_sqr_max = nqx * nqx * dqs[0] * dqs[0]
                    + nqy * nqy * dqs[1] * dqs[1]
                    + nqz * nqz * dqs[2] * dqs[2];
            real q_max = sqrt(q_sqr_max);
            int num_bins = (int) ceil(q_max / bin_size);
            vector<real> sq_bin;
            vector<real> q_bin;
            vector<int> count_bin;
            sq_bin.resize(num_bins);
            q_bin.resize(num_bins);
            count_bin.resize(num_bins);

            if (myrank == 0) {
                cout << nprocs << " CPUs\n\n"
                        << "bin size \t" << bin_size << "\n"
                        << "q_max    \t" << q_max << "\n";
            }

            real n_reci = 1. / num_part;
            real scos_local = 0; //will store cos-sum on each CPU
            real ssin_local = 0; //will store sin-sum on each CPU
            int ppp = (int) ceil((double) num_part / nprocs); //particles per proc

            Real3D coordP;

            python::list pyli;

            //loop over different q values
            //starting from zero because combinations with negative components 
            //will give the same result in S(q). so S(q) is the same for
            //the 8 vectors q=(x,y,z),(-x,y,z), (x,-y,z),(x,y,-z),(-x,-y,z),...
            for (int hx = -nqx; hx <= nqx; hx++) {
                for (int hy = -nqy; hy <= nqy; hy++) {
                    for (int hz = 0; hz <= nqz; hz++) {

                        //values of q-vector
                        q[0] = hx * dqs[0];
                        q[1] = hy * dqs[1];
                        q[2] = hz * dqs[2];
                        real q_abs = q.abs();

                        //determining the bin number
                        int bin_i = (int) floor(q_abs / bin_size);
                        q_bin[bin_i] += q_abs;
                        count_bin[bin_i] += 1;

                        //resetting the variables that store the local sum on each proc
                        scos_local = 0;
                        ssin_local = 0;

                        //loop over particles
                        for (int k = myrank * ppp; k < (1 + myrank) * ppp && k < num_part;
                                k++) {
                            coordP = config->getCoordinates(k);
                            scos_local += cos(q * coordP);
                            ssin_local += sin(q * coordP);
                        }
                        if (myrank != 0) {
                            boost::mpi::reduce(*system.comm, scos_local, plus<real > (), 0);
                            boost::mpi::reduce(*system.comm, ssin_local, plus<real > (), 0);
                        }

                        if (myrank == 0) {
                            real scos = 0;
                            real ssin = 0;
                            boost::mpi::reduce(*system.comm, scos_local, scos, plus<real > (), 0);
                            boost::mpi::reduce(*system.comm, ssin_local, ssin, plus<real > (), 0);
                            sq_bin[bin_i] += scos * scos + ssin * ssin;
                        }
                    }
                }
            }
            //creates the python list with the results            
            if (myrank == 0) {
                //starting with bin_i = 1 will leave out the value for q=0, otherwise start with bin_i=0
                for (int bin_i = 1; bin_i < num_bins; bin_i++) {
                    real c = (count_bin[bin_i]) ? 1 / (real) count_bin[bin_i] : 0;
                    sq_bin[bin_i] = n_reci * sq_bin[bin_i] * c;
                    q_bin[bin_i] = q_bin[bin_i] * c;

                    python::tuple q_Sq_pair;
                    q_Sq_pair = python::make_tuple(q_bin[bin_i], sq_bin[bin_i]);
                    pyli.append(q_Sq_pair);
                }
            }
            return pyli;
        }

        // this routine is for ordered configurations, e.g. particle 0 to 9 
        // belong to chain 1, particle 10 to 19 to chain 2 etc.        

        python::list StaticStructF::computeArraySingleChain(int nqx, int nqy, int nqz,
                real bin_factor, int chainlength) const {
            //fist the system coords are saved at each CPU
            System& system = getSystemRef();
            esutil::Error err(system.comm);
            Real3D Li = system.bc->getBoxL(); //Box size (Lx, Ly, Lz)

            int nprocs = system.comm->size(); // number of CPUs
            int myrank = system.comm->rank(); // current CPU's number

            int num_part = 0;
            ConfigurationPtr config = make_shared<Configuration > ();
            // loop over all CPU-numbers - to give all CPUs all particle coords
            for (int rank_i = 0; rank_i < nprocs; rank_i++) {
                map< size_t, Real3D > conf;
                if (rank_i == myrank) {
                    CellList realCells = system.storage->getRealCells();
                    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
                        int id = cit->id();
                        conf[id] = cit->position();
                    }
                }
                boost::mpi::broadcast(*system.comm, conf, rank_i);

                // for simplicity we will number the particles from 0
                for (map<size_t, Real3D>::iterator itr = conf.begin(); itr != conf.end(); ++itr) {
                    size_t id = itr->first;
                    Real3D p = itr->second;
                    config->set(id, p[0], p[1], p[2]);
                    //config->set(num_part, p[0], p[1], p[2]);
                    num_part++;
                }
            }
            cout << "particles are given to each CPU!\n";
            // now all CPUs have all particle coords and num_part is the total number
            // of particles   

            // use all CPUs
            // TODO it could be a problem if   n_nodes > num_part

            // here starts calculation of the static structure factor

            //step size for qx, qy, qz
            real dqs[3];
            dqs[0] = 2. * M_PIl / Li[0];
            dqs[1] = 2. * M_PIl / Li[1];
            dqs[2] = 2. * M_PIl / Li[2];

            Real3D q;

            //calculations for binning
            real bin_size = bin_factor * min(dqs[0], (dqs[1], dqs[2]));
            real q_sqr_max = nqx * nqx * dqs[0] * dqs[0]
                    + nqy * nqy * dqs[1] * dqs[1]
                    + nqz * nqz * dqs[2] * dqs[2];
            real q_max = sqrt(q_sqr_max);
            int num_bins = (int) ceil(q_max / bin_size);
            vector<real> sq_bin;
            vector<real> q_bin;
            vector<int> count_bin;
            sq_bin.resize(num_bins);
            q_bin.resize(num_bins);
            count_bin.resize(num_bins);

            if (myrank == 0) {
                cout << nprocs << " CPUs\n\n"
                        << "bin size \t" << bin_size << "\n"
                        << "q_max    \t" << q_max << "\n";
            }

            real n_reci = 1. / num_part;
            real chainlength_reci = 1. / chainlength;
            real scos_local = 0; //will store cos-sum on each CPU
            real ssin_local = 0; //will store sin-sum on each CPU   
            //will store the summation of the the single chain structure factor
            real singleChain_localSum = 0;
            Real3D coordP;
            python::list pyli;

            //calculations for parallelizing (over chains)
            int num_chains;
            if (num_part % chainlength == 0)
                num_chains = num_part / chainlength;
            else {
                cout << "ERROR: chainlenght does not match total number of "
                        << "particles. num_part % chainlenght is unequal 0. \n"
                        << "Calculation of SingleChain_StaticStructF aborted\n";
                return pyli;
            }
            int cpp = (int) ceil((double) num_chains / nprocs); //chains per proc
            cout << "chains per proc\t" << cpp << "\n";


            //loop over different q values
            //starting from zero because combinations with negative components 
            //will give the same result in S(q). so S(q) is the same for
            //the 8 vectors q=(x,y,z),(-x,y,z), (x,-y,z),(x,y,-z),(-x,-y,z),...
            for (int hx = -nqx; hx <= nqx; hx++) {
                for (int hy = -nqy; hy <= nqy; hy++) {
                    for (int hz = 0; hz <= nqz; hz++) {

                        //values of q-vector
                        q[0] = hx * dqs[0];
                        q[1] = hy * dqs[1];
                        q[2] = hz * dqs[2];
                        real q_abs = q.abs();

                        //determining the bin number
                        int bin_i = (int) floor(q_abs / bin_size);
                        q_bin[bin_i] += q_abs;
                        count_bin[bin_i] += 1;

                        //resetting the variable that stores the sum for each q-vector
                        singleChain_localSum = 0;

                        //loop over chains (cid is chain_id)
                        for (int cid = myrank * cpp; cid < (1 + myrank) * cpp
                                && cid < num_chains; cid++) {
                            scos_local = 0; //resetting the cos sum for the each chain
                            ssin_local = 0; //resetting the sin sum for the each chain
                            //loop over particles
                            for (int k = cid * chainlength; k < (1 + cid) * chainlength && k < num_part;
                                    k++) {
                                coordP = config->getCoordinates(k);
                                scos_local += cos(q * coordP);
                                ssin_local += sin(q * coordP);
                            }
                            //the (summation part of the) single chain structure 
                            // factors are summed up for the averaging at the 
                            // end (over the chains)
                            singleChain_localSum += scos_local * scos_local
                                    + ssin_local * ssin_local;
                        }


                        if (myrank != 0) {
                            boost::mpi::reduce(*system.comm, singleChain_localSum, plus<real > (), 0);
                        }

                        if (myrank == 0) {
                            real singleChainSum = 0;
                            boost::mpi::reduce(*system.comm, singleChain_localSum, singleChainSum, plus<real > (), 0);
                            sq_bin[bin_i] += singleChainSum;
                        }
                    }
                }
            }
            //creates the python list with the results            
            if (myrank == 0) {
                //starting with bin_i = 1 will leave out the value for q=0, otherwise start with bin_i=0
                for (int bin_i = 1; bin_i < num_bins; bin_i++) {
                    real c = (count_bin[bin_i]) ? 1 / (real) count_bin[bin_i] : 0;
                    sq_bin[bin_i] = n_reci * chainlength_reci * sq_bin[bin_i] * c;
                    q_bin[bin_i] = q_bin[bin_i] * c;

                    python::tuple q_Sq_pair;
                    q_Sq_pair = python::make_tuple(q_bin[bin_i], sq_bin[bin_i]);
                    pyli.append(q_Sq_pair);
                }
            }
            return pyli;
        }

        // TODO: this dummy routine is still needed as we have not yet ObservableVector
        // there has to be a function 'compute' because of the used template
        // otherwise a compiling error will occur

        real StaticStructF::compute() const {
            return -1.0;
        }

        void StaticStructF::registerPython() {
            using namespace espressopp::python;
            class_<StaticStructF, bases< Observable > >
                    ("analysis_StaticStructF", init< shared_ptr< System > >())
                    .def("compute", &StaticStructF::computeArray)
                    .def("computeSingleChain", &StaticStructF::computeArraySingleChain)
                    ;
        }
    }
}
