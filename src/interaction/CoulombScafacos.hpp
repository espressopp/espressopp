/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2022
      Data Center, Johannes Gutenberg University Mainz

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
#ifndef _INTERACTION_COULOMBSCAFACOS_HPP
#define _INTERACTION_COULOMBSCAFACOS_HPP

#include <cmath>

#include <boost/signals2.hpp>
#include "mpi.hpp"
#include "Potential.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"
#include "iterator/CellListIterator.hpp"
#include "Cell.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

#include "Tensor.hpp"

#include "System.hpp"

#include "iostream"
#include <fcs_config.h>
#include <cstdlib>
#include <cstring>
#include "fcs.h"
//#include "common/near/near.h"

#include "boost/serialization/vector.hpp"
#include "boost/serialization/complex.hpp"

using namespace std;

typedef complex<espressopp::real> dcomplex;

// the following two constants are not defined everywhere (e.g. not in Mac OS X)
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_2_SQRTPIl
#define M_2_SQRTPIl 1.1283791670955125738961589031215452L
#endif

#define M_2PI (2 * M_PIl)
#define M_PI2 (M_PIl * M_PIl)

#define M_1_SQRTPI (M_2_SQRTPIl * 0.5) /* 2/sqrt(pi)/2 = 1/sqrt(pi) */

extern int ifTuned;
extern fcs_int p3m_grid, p3m_cao;
extern fcs_float p3m_rcut, p3m_alpha;

namespace espressopp
{
namespace interaction
{
/** This class provides methods to compute forces and energies of the
 *  CoulombScafacos part. Currently it works with cubes and rectangular cuboids.
 *  Does not work for triclinic box, slab geometry.
 */

// TODO should be optimized (force energy and virial calculate the same stuff)

class CoulombScafacos : public PotentialTemplate<CoulombScafacos>
{
private:
    real prefactor;
    real tolerance;  // Ewald parameter
    int ntotal;      // cutoff in k space
    std::string method;

    shared_ptr<System> system;  // we need the system object to be able to access the box
                                // dimensions, communicator, number of particles, signals

    real Lx, Ly, Lz;  // local variable for system size
    int nParticles;   // local variable for the number of particles
    int num_glob;     // the total number of particles, input from arg/env
    bool ifBoundaryCross;

    // real rclx, rcly, rclz;
    real force_prefac[3];  // real array for force prefactors [0]: x,[1]: y,[2]: z

    // precalculated factors for the virial and virial tensor
    vector<real> virialPref;
    vector<Tensor> virialTensorPref;
    Tensor I;

    real sum_q2;
    real en_local;

    dcomplex *sum;
    dcomplex *totsum;

    // MPI_Comm communicator = MPI_COMM_WORLD;

    fcs_float *sfcs_coor;
    fcs_float *sfcs_cg;
    fcs_float *sfcs_field;
    fcs_float *sfcs_pot;
    fcs_float sfcs_prefac;
    fcs_int *sfcs_map;
    // fcs_int sfcs_pdc[3];
#if FCS_ENABLE_DIRECT
    char DIRECT_parameters[200];
#endif
#if FCS_ENABLE_EWALD
    char EWALD_parameters[200];
#endif
#if FCS_ENABLE_FMM
    char FMM_parameters[200];
#endif
#if FCS_ENABLE_MMM1D
    char MMM1D_parameters[200];
#endif
#if FCS_ENABLE_MMM2D
    char MMM2D_parameters[200];
#endif
#if FCS_ENABLE_PEPC
    char PEPC_parameters[200];
#endif
#if FCS_ENABLE_PP3MG
    char PP3MG_parameters[200];
#endif
#if FCS_ENABLE_P2NFFT
    char P2NFFT_parameters[200];
#endif
#if FCS_ENABLE_P3M
    char P3M_parameters[200];
#endif
#if FCS_ENABLE_VMG
    char VMG_parameters[200];
#endif

    fcs_int num_local;

    FCS handle;
    FCSResult result;

    char common_parameters[250];

public:
    static void registerPython();

    CoulombScafacos(shared_ptr<System> _system,
                    real _prefactor,
                    real _tolerance,
                    int _ntotal,
                    const char *_method);

    ~CoulombScafacos();

    // at this point we are ready to prepare the kvector[], it can be done just once at the begin
    void preset()
    {
        // TODO it could be parallelized too
        mpi::communicator communic = *system->comm;

        Real3D Li = system->bc->getBoxL();  // getting the system size
        Lx = Li[0];
        Ly = Li[1];
        Lz = Li[2];

        string Llx = to_string(Lx);
        string Lly = to_string(Ly);
        string Llz = to_string(Lz);
        string fcs_tol_field = to_string(tolerance);

        const char *llx = Llx.c_str();
        const char *lly = Lly.c_str();
        const char *llz = Llz.c_str();
        const char *ffcs_tol_field = fcs_tol_field.c_str();

        num_glob = ntotal;
        ifBoundaryCross = false;
        // MPI_Comm communicator = MPI_COMM_WORLD;
        // cout<<"SCATEST: 01 "<<num_glob<<"\n";
        sfcs_prefac = prefactor;  // 332.07;
                                  // fcs_int sfcs_pdc[3];
#if FCS_ENABLE_DIRECT
        strcpy(DIRECT_parameters, "direct_cutoff,0.0,direct_periodic_images,1,1,1");
#endif
#if FCS_ENABLE_EWALD
        strcpy(EWALD_parameters, "ewald_required_accuracy,1e-6");
#endif
#if FCS_ENABLE_FMM
        /*
        fcs_int FMM_absrel = FCS_METHOD_FMM_CUSTOM_RELATIVE;
        fcs_int FMM_dipole_correction = FCS_METHOD_FMM_STANDARD_DIPOLE_CORRECTION;
        fcs_float FMM_deltaE = 1e-6;
        */
        strcpy(FMM_parameters,
               "fmm_absrel,2,fmm_dipole_correction,0,fmm_tolerance_energy,1e-6,fmm_internal_tuning,"
               "0ll");
#endif
#if FCS_ENABLE_MMM1D
        strcpy(MMM1D_parameters,
               "mmm1d_bessel_cutoff,3,mmm1d_far_switch_radius,6.0,mmm1d_maxPWerror,1e-6");
#endif
#if FCS_ENABLE_MMM2D
        strcpy(MMM2D_parameters,
               "mmm2d_maxPWerror,1e-3,mmm2d_far_cutoff,1.73,mmm2d_dielectric_contrasts,3.17,2.13,"
               "mmm2d_layers_per_node,100,mmm2d_skin,0.5");
#endif
#if FCS_ENABLE_PEPC
        /*
        fcs_float PEPC_epsilon = 0.5;
        fcs_float PEPC_theta = 0.5;
        fcs_int PEPC_debuglevel = -1;
        */
        strcpy(PEPC_parameters, "pepc_debuglevel,-1,pepc_epsilon,0.5,pepc_theta,0.5");
#endif
#if FCS_ENABLE_PP3MG
        /*
        int *PP3MG_dims;
        fcs_int PP3MG_cells_x = 128;
        fcs_int PP3MG_cells_y = 128;
        fcs_int PP3MG_cells_z = 128;
        fcs_int PP3MG_ghost_cells = 4;
        fcs_int PP3MG_pol_degree = 4;
        fcs_float PP3MG_err_bound = 1e-6;
        fcs_int PP3MG_max_iterations = 20;
        fcs_int PP3MG_max_particles = 20;
        fcs_int PP3MG_periodic = 1;
        */
        strcpy(PP3MG_parameters,
               "pp3mg_cells_x,128,pp3mg_cells_y,128,pp3mg_cells_z,128,pp3mg_ghosts,4");
#endif
#if FCS_ENABLE_P2NFFT
        /* fcs_float P2NFFT_accuracy = 1e-6; */
        strcpy(P2NFFT_parameters, "p2nfft_required_accuracy,1e-6");
#endif
#if FCS_ENABLE_P3M
        /* fcs_float P3M_accuracy = 1e-6; */
        strcpy(P3M_parameters, "p3m_tolerance_field_abs,1e-6");
#endif
#if FCS_ENABLE_VMG
        strcpy(VMG_parameters,
               "vmg_cycle_type,2,vmg_max_iterations,20,vmg_max_level,6,vmg_near_field_cells,6,vmg_"
               "precision,1e-6,vmg_smoothing_steps,3");
#endif
        // if (handle == NULL)
        //  ;
        // else{
        //  fcs_destroy(handle);
        handle = NULL;
        //}
        result = NULL;

        // common_parameters="box_a,50.00,0.0,0.0,box_b,0.0,50.00,0.0,box_c,0.0,0.0,50.00,periodicity,1,1,1,offset,0.0,0.0,0.0,near_field_flag,1");
        strcpy(common_parameters, "box_a,");
        strcat(common_parameters, llx);
        strcat(common_parameters, ",0.0");
        strcat(common_parameters, ",0.0");
        strcat(common_parameters, ",box_b");
        strcat(common_parameters, ",0.0,");
        strcat(common_parameters, lly);
        strcat(common_parameters, ",0.0");
        strcat(common_parameters, ",box_c");
        strcat(common_parameters, ",0.0");
        strcat(common_parameters, ",0.0,");
        strcat(common_parameters, llz);

        strcat(common_parameters, ",periodicity,1,1,1");
        strcat(common_parameters, ",offset,0.0,0.0,0.0");
        strcat(common_parameters, ",near_field_flag,1");
        strcat(common_parameters, ",tolerance_field,");
        strcat(common_parameters, ffcs_tol_field);

        // getParticleNumber();
        // num_local=nParticles;
        // //MPI_Comm_size (communicator, &comm_size);
        // //fprintf (stderr, "comm_size: %d \n", comm_size);
        // //MPI_Dims_create (comm_size, 3, dims);
        // //fprintf (stderr, "dims: %d %d %d\n", dims[0], dims[1], dims[2]);
        // //MPI_Cart_create (communicator, 3, dims, pers, 0, &cart_comm);
        // //communicator = cart_comm;
        // sfcs_coor = (fcs_float *) malloc (3 * sizeof (fcs_float) * num_local);
        // sfcs_cg = (fcs_float *) malloc (sizeof (fcs_float) * num_local);
        // sfcs_field = (fcs_float *) malloc (3 * sizeof (fcs_float) * num_local);
        // sfcs_pot = (fcs_float *) malloc (sizeof (fcs_float) * num_local);
        // //velocities = (fcs_float *) malloc (3 * sizeof (fcs_float) * num_local);
        // //masses = (fcs_float *) malloc (sizeof (fcs_float) * num_local);

        result = fcs_init(&handle, method.c_str(), communic);
        // fcs_result_print_result (result);
        fcs_result_destroy(result);
        communic.barrier();

        result = fcs_set_parameters(handle, common_parameters, FCS_FALSE);
        // fcs_result_print_result (result);
        fcs_result_destroy(result);

        result = fcs_set_total_particles(handle, num_glob);  // NOT LOCAL
        // fcs_result_print_result (result);
        fcs_result_destroy(result);
        communic.barrier();

        switch (fcs_get_method(handle))
        {
#if FCS_ENABLE_DIRECT
            case FCS_METHOD_DIRECT:
                result = fcs_set_parameters(handle, DIRECT_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_EWALD
            case FCS_METHOD_EWALD:
                result = fcs_set_parameters(handle, EWALD_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_FMM
            case FCS_METHOD_FMM:
                result = fcs_set_parameters(handle, FMM_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_MMM1D
            case FCS_METHOD_MMM1D:
                result = fcs_set_parameters(handle, MMM1D_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_MMM2D
            case FCS_METHOD_MMM2D:
                result = fcs_set_parameters(handle, MMM2D_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_PEPC
            case FCS_METHOD_PEPC:
                result = fcs_set_parameters(handle, PEPC_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_PP3MG
            case FCS_METHOD_PP3MG:
                result = fcs_set_parameters(handle, PP3MG_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_P2NFFT
            case FCS_METHOD_P2NFFT:
                result = fcs_set_parameters(handle, P2NFFT_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif
#if FCS_ENABLE_P3M
            case FCS_METHOD_P3M:
                result = fcs_set_parameters(handle, P3M_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                ifBoundaryCross = true;
                break;
#endif
#if FCS_ENABLE_VMG
            case FCS_METHOD_VMG:
                result = fcs_set_parameters(handle, VMG_parameters, FCS_FALSE);
                fcs_result_destroy(result);
                break;
#endif

            default:
                // result = fcs_set_parameters( handle, EWALD_parameters, FCS_FALSE);
                // fcs_result_destroy (result);
                break;
        }
        communic.barrier();

        if (sum != NULL)
        {
            delete[] sum;
            sum = NULL;
        }
        if (totsum != NULL)
        {
            delete[] totsum;
            totsum = NULL;
        }
    }

    // here we get the current particle number on the current node
    // and set the auxiliary arrays eikx, eiky, eikz
    void getParticleNumber() { nParticles = system->storage->getNRealParticles(); }

    // it counts the squared charges over all system. It is used for self energy calculations
    void count_charges(CellList realcells)
    {
        real node_sum_q2 = 0.0;
        for (iterator::CellListIterator it(realcells); !it.isDone(); ++it)
        {
            Particle &p = *it;
            node_sum_q2 += (p.q() * p.q());
        }
        sum_q2 = 0.0;
        mpi::all_reduce(*system->comm, node_sum_q2, sum_q2, plus<real>());
    }

    // set/get the parameters
    void setPrefactor(real _prefactor)
    {
        prefactor = _prefactor;
        preset();
    }
    real getPrefactor() const { return prefactor; }
    void setTolerance(real _tolerance)
    {
        tolerance = _tolerance;
        preset();
    }
    real getTolerance() const { return tolerance; }
    void setNTotal(int _ntotal)
    {
        ntotal = _ntotal;
        preset();
    }
    int getNTotal() const { return ntotal; }
    void setMethod(const char *_method)
    {
        method = _method;
        preset();
    }
    char const *getMethod() const { return method.c_str(); }

    real _computeEnergy(CellList realcells)
    {
        mpi::communicator communic = *system->comm;

        // int n_nodes = communic.size();
        // int this_node = communic.rank();

        // real fact;
        real energy = 0;

        mpi::all_reduce(communic, en_local, energy, plus<real>());

        return energy;
    }

    // @TODO this function could be void,
    bool _computeForce(CellList realcells)
    {
        mpi::communicator communic = *system->comm;

        // fcs_int i, j;
        fcs_int k, k3, idx, nlocal_map;

        getParticleNumber();
        num_local = nParticles;
        // cout<<"SCATEST: NUMLC = "<<num_local<<" "<<Lx<<" "<<Ly<<" "<<Lz<<endl;
        // MPI_Comm_size (communicator, &comm_size);
        // fprintf (stderr, "comm_size: %d \n", comm_size);
        // MPI_Dims_create (comm_size, 3, dims);
        // fprintf (stderr, "dims: %d %d %d\n", dims[0], dims[1], dims[2]);
        // MPI_Cart_create (communicator, 3, dims, pers, 0, &cart_comm);
        // communicator = cart_comm;
        sfcs_coor = (fcs_float *)malloc(3 * sizeof(fcs_float) * num_local);
        sfcs_cg = (fcs_float *)malloc(sizeof(fcs_float) * num_local);
        sfcs_field = (fcs_float *)malloc(3 * sizeof(fcs_float) * num_local);
        sfcs_pot = (fcs_float *)malloc(sizeof(fcs_float) * num_local);
        // velocities = (fcs_float *) malloc (3 * sizeof (fcs_float) * num_local);
        // masses = (fcs_float *) malloc (sizeof (fcs_float) * num_local);
        sfcs_map = (fcs_int *)malloc(sizeof(fcs_int) * num_local);
        // cout<<"SCATEST: DEF end\n";
        // load all particles from the current processor and transform into
        // fcs format
        k = 0;
        idx = 0;
        if (ifBoundaryCross)
        {
            for (iterator::CellListIterator it(realcells); !it.isDone(); ++it)
            {
                Particle &p = *it;

                if (p.q() != 0.0)
                {
                    k3 = k * 3;

                    sfcs_coor[k3] = p.position()[0];
                    sfcs_coor[k3 + 1] = p.position()[1];
                    sfcs_coor[k3 + 2] = p.position()[2];

                    if (sfcs_coor[k3] < 0.0)
                        sfcs_coor[k3] += Lx;
                    else if (sfcs_coor[k3] > Lx)
                        sfcs_coor[k3] -= Lx;
                    if (sfcs_coor[k3 + 1] < 0.0)
                        sfcs_coor[k3 + 1] += Ly;
                    else if (sfcs_coor[k3 + 1] > Ly)
                        sfcs_coor[k3 + 1] -= Ly;
                    if (sfcs_coor[k3 + 2] < 0.0)
                        sfcs_coor[k3 + 2] += Lz;
                    else if (sfcs_coor[k3 + 2] > Lz)
                        sfcs_coor[k3 + 2] -= Lz;
                    sfcs_cg[k] = p.q();

                    sfcs_map[k] = idx;
                    k++;
                }

                idx++;
            }

            nlocal_map = k;
        }
        else
        {
            for (iterator::CellListIterator it(realcells); !it.isDone(); ++it)
            {
                Particle &p = *it;

                if (p.q() != 0.0)
                {
                    k3 = k * 3;
                    sfcs_coor[k3] = p.position()[0];
                    sfcs_coor[k3 + 1] = p.position()[1];
                    sfcs_coor[k3 + 2] = p.position()[2];
                    sfcs_cg[k] = p.q();

                    sfcs_map[k] = idx;
                    k++;
                }

                idx++;
            }

            nlocal_map = k;
        }
        if (ifTuned)
        {
            fcs_p3m_get_alpha(handle, &p3m_alpha);
            fcs_p3m_get_r_cut(handle, &p3m_rcut);
            fcs_p3m_get_cao(handle, &p3m_cao);
            fcs_p3m_get_grid(handle, &p3m_grid);
            // fcs_p3m_get_total_energy(handle,&tot_en);
            // result=fcs_tune(handle, num_local, sfcs_coor, sfcs_cg);
        }
        else
        {
            result = fcs_tune(handle, nlocal_map, sfcs_coor, sfcs_cg);
            fcs_p3m_get_alpha(handle, &p3m_alpha);
            fcs_p3m_get_r_cut(handle, &p3m_rcut);
            fcs_p3m_get_cao(handle, &p3m_cao);
            fcs_p3m_get_grid(handle, &p3m_grid);
            // fcs_p3m_get_total_energy(handle,&tot_en);
            ifTuned = 1;
        }

        // fcs_result_print_result (result);
        // cout<<"Destroy... \n";
        fcs_result_destroy(result);

        // cout<<"Barrier... \n";
        communic.barrier();
        // cout<<" TUN_RES: \n";
        // fcs_result_print_result (result);
        // cout<<" PRM: \n";
        // fcs_print_parameters (handle);

        // cout<<" RUN: \n";
        result = fcs_run(handle, nlocal_map, sfcs_coor, sfcs_cg, sfcs_field, sfcs_pot);
        communic.barrier();

        // fcs_destroy (handle);

        // cout<<" CR: "<<sfcs_coor[0]<<" "<<sfcs_coor[1]<<" "<<sfcs_coor[2]<<"\n";
        // cout<<" FF: "<<sfcs_field[0]<<" "<<sfcs_field[1]<<" "<<sfcs_field[2]<<"\n";
        // cout<<" PREF= "<<sfcs_prefac<<endl;
        k = 0, idx = 0;

        for (iterator::CellListIterator it(realcells); !it.isDone(); ++it)
        {
            Particle &p = *it;

            if (sfcs_map[k] == idx)
            {
                k3 = k * 3;
                p.force()[0] += sfcs_prefac * sfcs_field[k3] * p.q();
                p.force()[1] += sfcs_prefac * sfcs_field[k3 + 1] * p.q();
                p.force()[2] += sfcs_prefac * sfcs_field[k3 + 2] * p.q();

                k++;
            }

            idx++;
        }

        for (k = 0, en_local = .0; k < nlocal_map; k++)
        {
            en_local += sfcs_prefac * sfcs_pot[k] * sfcs_cg[k] / 2.0;
        }

        free(sfcs_coor);
        free(sfcs_cg);
        free(sfcs_field);
        free(sfcs_pot);
        free(sfcs_map);

        return true;
    }

    // compute virial for this interaction
    // (!note: all particle interaction contains only one potential)
    real _computeVirial(CellList realcells)
    {
        mpi::communicator communic = *system->comm;

        // int n_nodes = communic.size();
        // int this_node = communic.rank();

        real node_virial = 0;
        real virial = 0;
        mpi::all_reduce(communic, node_virial, virial, plus<real>());

        return virial;
    }

    // compute virial Tensor for this interaction
    // (!note: all particle interaction contains only one potential)
    Tensor _computeVirialTensor(CellList realcells)
    {
        mpi::communicator communic = *system->comm;

        // int n_nodes = communic.size();
        // int this_node = communic.rank();

        Tensor node_virialTensor = 0;

        Tensor virialTensor(0.0);
        mpi::all_reduce(communic, node_virialTensor, virialTensor, plus<Tensor>());

        // using boost::mpi::all_reduce or reduce is very slow for <Tensor>
        // as a suggestion, one could try the line below:
        // mpi::all_reduce( communic, (double*)&node_virialTensor, 6, (double*)&virialTensor,
        // plus<double>());

        return virialTensor;
    }

    real _computeEnergySqrRaw(real distSqr) const
    {
        esutil::Error err(system->comm);
        stringstream msg;
        msg << "There is no sense to call this function for Ewald summation";
        err.setException(msg.str());
        return 0.0;
    }
    bool _computeForceRaw(Real3D &force, const Real3D &dist, real distSqr) const
    {
        esutil::Error err(system->comm);
        stringstream msg;
        msg << "There is no sense to call this function for Ewald summation";
        err.setException(msg.str());
        return false;
    }
};
}  // namespace interaction
}  // namespace espressopp

#endif
