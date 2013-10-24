#include "python.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "Configuration.hpp"
#include "StaticStructF.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"

#include <boost/serialization/map.hpp>

#include <math.h>       //cos and ceil and sqrt
#include <algorithm>    // std::min
#include <functional>   // std::plus
#include <fstream> //outputfile
#include <iomanip> //setprecision

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

using namespace espresso;
using namespace espresso::iterator;
using namespace std;

namespace espresso {
    namespace analysis {

        // TODO currently works correctly for Lx = Ly = Lz
        // nqx is a number which corresponds to the different x-values of the
        // diffraction vector q. greater nqx produces more different x-values   

        python::list StaticStructF::computeArray(int nqx, int nqy, int nqz,
                real bin_factor) const {
            //fist the system coords are saved at each CPU
            System& system = getSystemRef();
            esutil::Error err(system.comm);
            Real3D Li = system.bc->getBoxL(); //Box size (Lx, Ly, Lz)
            //Real3D Li_half = Li / 2.; //FM dont need for staticstructf

            int nprocs = system.comm->size(); // number of CPUs
            int myrank = system.comm->rank(); // current CPU's number
            //cout << "1 \n"; //flag
            // FM do not need this - no histogram
            /*
            real histogram [rdfN];
            for(int i=0;i<rdfN;i++) histogram[i]=0;
        
            real dr = Li_half[1] / (real)rdfN; // If you work with nonuniform Lx, Ly, Lz, you
            // should use for Li_half[XXX] the shortest side length   
             */
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
                    //cout << "2 \n";   //flag 
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
            //cout << "3 \n";//flag
            // now all CPUs have all particle coords and num_part is the total number
            // of particles

            // use all cpus
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
            real scos = 0;
            real ssin = 0;
            real scos_local = 0;
            real ssin_local = 0;
            int ppp = (int) ceil((double) num_part / nprocs);
            //cout << "ppp " << ppp << "\n";
            Real3D coordP;

            python::list pyli;


            //loop over different q values
            for (int hx = -nqx; hx <= nqx; hx++) {
                for (int hy = -nqy; hy <= nqy; hy++) {
                    for (int hz = -nqz; hz <= nqz; hz++) {
                        //int hy = 1;
                        //int hz = 1;

                        //values of q-vector
                        q[0] = hx * dqs[0];
                        q[1] = hy * dqs[1];
                        q[2] = hz * dqs[2];
                        real q_abs = q.abs();
                        if (myrank == 0) {
                            cout << "q.abs for (" << hx << "," << hy << "," << hz << "): "
                                    << q_abs << "\n";
                        }
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

                        boost::mpi::reduce(*system.comm, scos_local, scos, plus<real > (), 0);
                        boost::mpi::reduce(*system.comm, ssin_local, ssin, plus<real > (), 0);

                        if (myrank == 0) {
                            sq_bin[bin_i] += scos * scos + ssin * ssin;
                        }
                    }
                }
            }
            //creates an output file with q and S(q) values
            ofstream outfile;
            if (myrank == 0) {
                outfile.open("q_sq_values.txt");
                for (int bin_i = 1; bin_i < num_bins; bin_i++) {
                    real c = (count_bin[bin_i]) ? 1 / (real) count_bin[bin_i] : 0;
                    sq_bin[bin_i] = n_reci * sq_bin[bin_i] * c;
                    q_bin[bin_i] = q_bin[bin_i] * c;
                    pyli.append(sq_bin[bin_i]);
                    
                    if(outfile.is_open()){
                        outfile << setprecision(8);
                        outfile << fixed;
                        outfile << q_bin[bin_i] << "   \t" << sq_bin[bin_i] << "\n";                       
                    }
                    else cout << "Unable to open output file";                    
                }
                outfile.close();
            }

            return pyli;

        }
        // TODO: this dummy routine is still needed as we have not yet ObservableVector

        real StaticStructF::compute() const {
            return -1.0;
        }

        void StaticStructF::registerPython() {
            using namespace espresso::python;
            class_<StaticStructF, bases< Observable > >
                    ("analysis_StaticStructF", init< shared_ptr< System > >())
                    .def("compute", &StaticStructF::computeArray)
                    ;
        }
    }
}
