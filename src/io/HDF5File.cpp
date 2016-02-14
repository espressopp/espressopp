/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & Johannes Gutenberg-Universität Mainz

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


// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "HDF5File.hpp"  // keep python.hpp on top
#include "storage/Storage.hpp"
#include "System.hpp"
#include "storage/DomainDecomposition.hpp"
#include "bc/BC.hpp"
#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"
#include "iterator/CellListIterator.hpp"

#include "boost/mpi.hpp"
#include "boost/mpi/communicator.hpp"
#include "mpi.h"

#include <fstream>
#include <sstream>

#ifdef HDF5_LAYER
    #include "hdf5.h"
    #include "hdf5_hl.h"
#endif
//#include "H5Cpp.h"

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;


namespace espressopp {
  namespace io {

  	//using namespace iterator;

    void HDF5File::write_n_to_1(){


    	shared_ptr<System> system = getSystem();

    	char *ch_f_name = new char[file_name.length() + 1];
    	strcpy(ch_f_name, file_name.c_str());

    	int rank = system->comm->rank();
    	int mpi_ranks = system->comm->size();
    	string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();

    	//if (system->CommunicatorIsInitialized)

    	int ierr;

		int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
		int maxN;   // maximal number of particles one processor has
		int totalN; // total number of particles all processors have

		//std::cout << "part_rank[" << rank << "] = " << myN << std::endl;

		boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
		boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace

		int* array_nparticles = new int [mpi_ranks];   // to write contiguos in file

		boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguos writing

		/*
		for(int k = 0; k< mpi_ranks; k++) {


			std::cout << "part[" << k << "] = " << array_nparticles[k] << ",";
		}

		std::cout << std::endl;
		*/

		//cout << "rank: " << rank << " has " << maxN << " max particles" << endl;
		//cout << "rank: " << rank << " has " << myN << " particles" << endl;

		char dataSetName[256];
		int RANK = 2;


		typedef struct {
			size_t pid;
			size_t type;
			real mass;
			real charge;
			real lambda; // adress flag: 0 means no adress
			real drift;
			real lambdaDeriv;
			int state;
			double x[3];
			double v[3];
			double force[3];
		} particle_info;

		CellList realCells = system->storage->getRealCells();
		//MPI_Info info;
		MPI_Info info  = MPI_INFO_NULL;
		//MPI_Info_create(&info);
		//MPI_Info_set(info, "file", "jack");

		// create type for array-like objects, like coordinates, vel and force
		hsize_t dimearr[1] = {3};
		hid_t loctype = H5Tarray_create1(H5T_NATIVE_DOUBLE, 1, dimearr, NULL);


		//std::cout << "Create HDF5 type" << std::endl;
		/*
		hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(record_t));
		H5Tinsert(particle_record, "pid", HOFFSET(record_t,pid), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "type", HOFFSET(record_t,type), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "x", HOFFSET(record_t,x), loctype);
		H5Tinsert(particle_record, "v", HOFFSET(record_t,v), loctype);
		*/

		// create the HDF5 compound datatype from the struct
		hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(particle_info));
		H5Tinsert(particle_record, "pid", HOFFSET(particle_info,pid), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "type", HOFFSET(particle_info,type), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "mass", HOFFSET(particle_info,mass), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "charge", HOFFSET(particle_info,charge), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "lambda", HOFFSET(particle_info,lambda), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "drift", HOFFSET(particle_info,drift), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "lambdaDeriv", HOFFSET(particle_info,lambdaDeriv), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "state", HOFFSET(particle_info,state), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "x", HOFFSET(particle_info,x), loctype);
		H5Tinsert(particle_record, "v", HOFFSET(particle_info,v), loctype);
		H5Tinsert(particle_record, "force", HOFFSET(particle_info,force), loctype);



		hid_t acc_template;
		hid_t file_id, dset_id;
		hid_t filespace, memspace;
		herr_t status;

		hsize_t count[RANK];
		hsize_t offset[RANK];

		hsize_t dimsf[RANK];


		acc_template = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);
		file_id = H5Fcreate(ch_f_name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
		assert(file_id > 0);

		H5Pclose(acc_template);

		//dataspace = H5Screate_simple(RANK, dimsf, NULL);



		/*
		 * assert() section before filling the data
		 *
		 */

		//assert( maxN == myN);

		//assert(offset[0] == (rank * myN));
		//assert(offset[1] == 0);

		particle_info* particles_u  = new particle_info [myN];

		int i = 0;
		assert( i == 0);

		  if( unfolded ){
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();
			  Int3D& img = cit->image();
			  Real3D L = system->bc->getBoxL();
			  particles_u[i].pid = cit->id();
			  particles_u[i].type = cit->type();
			  particles_u[i].mass = cit->mass();
			  particles_u[i].charge = cit->q();
			  particles_u[i].lambda = cit->lambda();
			  particles_u[i].drift = cit->drift();
			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
			  particles_u[i].state = cit->state();
			  particles_u[i].x[0] = pos[0] + img[0] * L[0];
			  particles_u[i].x[1] = pos[1] + img[1] * L[1];
			  particles_u[i].x[2] = pos[2] + img[2] * L[2];
			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];
			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];

			  i++;
			}
		  }
		  else{
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();
			  particles_u[i].pid = cit->id();
			  particles_u[i].type = cit->type();
			  particles_u[i].mass = cit->mass();
			  particles_u[i].charge = cit->q();
			  particles_u[i].lambda = cit->lambda();
			  particles_u[i].drift = cit->drift();
			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
			  particles_u[i].state = cit->state();
			  particles_u[i].x[0] = pos[0];
			  particles_u[i].x[1] = pos[1];
			  particles_u[i].x[2] = pos[2];
			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];
			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];

			  i++;
			}
		  }
		 // std::cout << "Finished filling the particle records!" << std::endl;

        hsize_t di[RANK];
        //di[0] = maxN*8;
        di[0] = myN;
        di[1] = 6;

        hsize_t di2[RANK];
		di2[0] = myN;
		di2[1] = H5Tget_size(particle_record);


        /** prepare chunking according to the number of particles myN in a rank **/

        hsize_t chunks[RANK];
        chunks[0] = myN;
        chunks[1] = 6;

        dimsf[0] = totalN;
        dimsf[1] = 1;


        /* Create a dataset creation property list and set it to use chunking
             * with a chunk size of 10x1 */
        //hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        //H5Pset_chunk(dcpl_id, RANK, chunks);

		filespace = H5Screate_simple(RANK, dimsf, NULL);

		//H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_dataspace, NULL, count_dataspace, NULL);

		//H5Sclose(dataspace);

		//memspace = H5Screate_simple(RANK, count_memspace, NULL);

        //H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_memspace, NULL, count_memspace, NULL);

		//std::string dataset_name = "PosVel_" + rankstring;
		//sprintf(dataSetName, dataset_name.c_str());
		sprintf(dataSetName, "Particles");
		dset_id = H5Dcreate2(file_id, dataSetName, particle_record, filespace,
																H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(filespace);

		int sumup = 0;
		count[0] = myN;
		count[1] = dimsf[1];
		if (rank == 0) {
		offset[0] = rank;
		} else {

			for(int L=0; L<rank; L++) {
				sumup += array_nparticles[L];
			}

			offset[0] = sumup;
		}

		offset[1] = 0;

		memspace = H5Screate_simple(RANK, count, NULL);


		filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);



		hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // to write datasets independently

		status = H5Dwrite(dset_id, particle_record, memspace, filespace, plist_id, particles_u);

		assert( status != -1);

		// check weather is the case to create scalar datasets to hold attribute!
		/*
		hid_t dataspace_id, attribute_id;
		hsize_t     dims;

		dims = 1;
		dataspace_id = H5Screate_simple(1, &dims, NULL);
		double timestep = integrator->getTimeStep();
		long long step = integrator->getStep();
		double attr_data_timestep = timestep;
		long long attr_data_step = step;

		attribute_id = H5Acreate(dset_id, "timestep", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT);

		status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data_timestep);
		attribute_id = H5Acreate(dset_id, "step", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data_step);
		H5Aclose(attribute_id);
		*/



		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Tclose(particle_record);
		H5Sclose(memspace);
		//H5Pclose(dcpl_id);


		//		status = H5Gclose (group);  //not generating this group terminates ok!!
		//			//H5Fflush(file);
	//}
		H5Pclose(plist_id);
		H5Fclose(file_id);
		//MPI_Info_free(&info);

		delete [] ch_f_name;

    	/** Subset write implementation? Probably not **/

		//free(coord);
		//delete [] b;
		//delete [] co;
		delete [] array_nparticles;
		delete [] particles_u;
	}




    void HDF5File::write_n_to_n(){

    	shared_ptr<System> system = getSystem();
    	//shared_ptr<integrator::MDIntegrator> integratori = integrator::;
    	//integrator->
		//ExtensionType exttype = integrator->getType();
    	//Extension(integrator);
    	//std::cout << "Entering Write routine" << std::endl;

    	int rank = system->comm->rank();

    	size_t filename_length = file_name.length();
    	string suffix = file_name.substr(filename_length-3, 3);

    	string base_filename = file_name.substr(0,filename_length-3);

    	//cout << "base filename: " << base_filename << std::endl;
    	//cout << "suffix: " << suffix << std::endl;

    	string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();

    	std::string final_name = base_filename + "_" + rankstring + suffix;

    	//char *ch_f_name = new char[file_name.length() + 1];
    	//strcpy(ch_f_name, file_name.c_str());
    	char *ch_f_name = new char[final_name.length() + 1];
    	strcpy(ch_f_name, final_name.c_str());

    	//if (system->CommunicatorIsInitialized)

    	/** N-to-N implementation **/

    	hid_t acc_template;
		hid_t dataspace, memspace, dset, file_id;
		herr_t status;
		int ierr;

		int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
		int maxN;   // maximal number of particles one processor has
		int totalN; // total number of particles all processors have


		char dataSetName[256];
		int RANK = 2;

		hsize_t start_dataspace[RANK];
		hsize_t stride_dataspace[RANK];
		hsize_t count_dataspace[RANK];
		hsize_t block_dataspace[RANK];

		start_dataspace[0] = rank*myN;
		start_dataspace[1] = 0;
		count_dataspace[0] = myN;
		count_dataspace[1] = 6;
		block_dataspace[0] = myN;
		block_dataspace[1] = 6;


		hsize_t start_memspace[RANK];
		hsize_t stride_memspace[RANK];
		hsize_t count_memspace[RANK];

		start_memspace[0] = 0;
		start_memspace[1] = 0;
		count_memspace[0] = myN;
		count_memspace[1] = 6;


		//double* coordinates = new real [3 * myN];  // buffer for gather
		//double* velocities  = new real [3 * myN];  // buffer for gather
		//int*  ids         = new int [myN];  // buffer for gather


		//double coord[count_memspace[0]][count_memspace[1]];

		CellList realCells = system->storage->getRealCells();

		// Define important information and build up compund datatype

		typedef struct {
		     size_t pid;   // pid
		     size_t type;   // type
		     double x[3]; // coordinates
		     double v[3];  // velocities

		  } record_t;


		typedef struct {
			size_t pid;
			size_t type;
			real mass;
			real charge;
			real lambda; // adress flag: 0 means no adress
			real drift;
			real lambdaDeriv;
			int state;
			double x[3];
			double v[3];
			double force[3];
		} particle_info;


		//typedef struct {
		//	double timo;
		//} timestep_info;



		  // 1 record_t represent a "particle"

		  hsize_t     dimearr[1] = {3};
		  hid_t loctype = H5Tarray_create1(H5T_NATIVE_DOUBLE, 1, dimearr, NULL);


		  //std::cout << "Create HDF5 type" << std::endl;
		  /*
		  hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(record_t));
		  H5Tinsert(particle_record, "pid", HOFFSET(record_t,pid), H5T_NATIVE_INT);
		  H5Tinsert(particle_record, "type", HOFFSET(record_t,type), H5T_NATIVE_INT);
		  H5Tinsert(particle_record, "x", HOFFSET(record_t,x), loctype);
		  H5Tinsert(particle_record, "v", HOFFSET(record_t,v), loctype);
		  */

		  // with extra particle info
		  hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(particle_info));
		  H5Tinsert(particle_record, "pid", HOFFSET(particle_info,pid), H5T_NATIVE_INT);
		  H5Tinsert(particle_record, "type", HOFFSET(particle_info,type), H5T_NATIVE_INT);
		  H5Tinsert(particle_record, "mass", HOFFSET(particle_info,mass), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "charge", HOFFSET(particle_info,charge), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "lambda", HOFFSET(particle_info,lambda), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "drift", HOFFSET(particle_info,drift), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "lambdaDeriv", HOFFSET(particle_info,lambdaDeriv), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "state", HOFFSET(particle_info,state), H5T_NATIVE_INT);

		  H5Tinsert(particle_record, "x", HOFFSET(particle_info,x), loctype);
		  H5Tinsert(particle_record, "v", HOFFSET(particle_info,v), loctype);
		  H5Tinsert(particle_record, "force", HOFFSET(particle_info,force), loctype);


		  //std::cout << "Created HDF5 type" << std::endl;

		  //std::cout << "Creating HDF5 file" << std::endl;

		file_id = H5Fcreate(ch_f_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		assert(file_id > 0);
		//std::cout << "Created HDF5 file" << std::endl;


		//std::cout << "Just before filling the particle records!" << std::endl;
		//record_t* particles_u  = new record_t [myN];
		particle_info* particles_u  = new particle_info [myN];

		int i = 0;
		assert( i == 0);

		  if( unfolded ){
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

			  //ids[i] = cit->id();

			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();
			  Int3D& img = cit->image();
			  Real3D L = system->bc->getBoxL();

			  particles_u[i].pid = cit->id();

			  particles_u[i].type = cit->type();

			  particles_u[i].mass = cit->mass();

			  particles_u[i].charge = cit->q();

			  particles_u[i].lambda = cit->lambda();

			  particles_u[i].drift = cit->drift();

			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();

			  particles_u[i].state = cit->state();

			  particles_u[i].x[0] = pos[0] + img[0] * L[0];
			  particles_u[i].x[1] = pos[1] + img[1] * L[1];
			  particles_u[i].x[2] = pos[2] + img[2] * L[2];

			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];

			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];



			  i++;
			}
		  }
		  else{
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

			  //ids[i] = cit->id();
			  //size_t tipo = cit->type();
			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();


			  particles_u[i].pid = cit->id();

			  particles_u[i].type = cit->type();

			  particles_u[i].mass = cit->mass();

			  particles_u[i].charge = cit->q();

			  particles_u[i].lambda = cit->lambda();

			  particles_u[i].drift = cit->drift();

			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();

			  particles_u[i].state = cit->state();

			  particles_u[i].x[0] = pos[0];
			  particles_u[i].x[1] = pos[1];
			  particles_u[i].x[2] = pos[2];

			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];

			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];


			  i++;
			}
		  }
		 // std::cout << "Finished filling the particle records!" << std::endl;


        hsize_t di[RANK];
        di[0] = myN;
        di[1] = 6;

        // new dimensions of dataspace
        hsize_t di2[RANK];
		di2[0] = myN;
		di2[1] = H5Tget_size(particle_record);


		//dataspace = H5Screate_simple(RANK, di, NULL);
		//new dataspace almost complete
		//dataspace = H5Screate_simple(RANK, di2, NULL);

		hsize_t test[1];
		test[0] = myN;

		dataspace = H5Screate_simple(1, test, NULL);
		//std::cout << "Created dataspace" << std::endl;



		sprintf(dataSetName, "Particles");
		//dset = H5Dcreate2(file_id, dataSetName, H5T_IEEE_F64BE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		// new dataset creation
		dset = H5Dcreate2(file_id, dataSetName, particle_record, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//std::cout << "Created dataset" << std::endl;

		//hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
		//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		//status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, coordinates);
		//status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, co);
		//status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, coord);
		//status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace, plist_id, coord);

		//status = H5Dwrite(dset, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
		//std::cout << "Before writing dataset" << std::endl;
		// new writes to dataset
		status = H5Dwrite(dset, particle_record, H5S_ALL, H5S_ALL, H5P_DEFAULT, particles_u);


		assert( status != -1);

		// create attribute for timestep/step

		/* Create the data space for the attribute. */

		hid_t dataspace_id, attribute_id;
		hsize_t     dims;

		dims = 1;
		dataspace_id = H5Screate_simple(1, &dims, NULL);
		double timestep = integrator->getTimeStep();
		long long step = integrator->getStep();
		double attr_data_timestep = timestep;
		long long attr_data_step = step;


	   /* Create a dataset attribute. */
		attribute_id = H5Acreate(dset, "timestep", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

	   /* Write the attribute data. */
		status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data_timestep);
		attribute_id = H5Acreate(dset, "step", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data_step);
		H5Aclose(attribute_id);
		H5Sclose(dataspace);
		//H5Sclose(memspace);
		H5Dclose(dset);

		// close type
		H5Tclose(particle_record);


		//		status = H5Gclose (group);  //not generating this group terminates ok!!
		//			//H5Fflush(file);
	//}
		//H5Pclose(plist_id);
		H5Fclose(file_id);
		//MPI_Info_free(&info);

		delete [] ch_f_name;

    	/** Subset write implementation? Probably not **/

		//free(coord);
		//delete [] b;
		//delete [] co;
		//delete [] coordinates;
		//delete [] velocities;
		delete [] particles_u;
		//delete [] ids;

      }



    void HDF5File::write(){

    	int iomodus = getIomode();
    	if (iomodus == 1 || iomodus == 0) {
    		write_n_to_1();
		} else if (iomodus == 2) {
			write_n_to_n();
		}
    }


    // Python wrapping
    void HDF5File::registerPython() {

      using namespace espressopp::python;

      class_<HDF5File, bases<ParticleAccess>, boost::noncopyable >
      ("io_HDF5File", init< shared_ptr< System >,
                           shared_ptr< integrator::MDIntegrator >,
                           std::string,
						   int,
                           bool,
                           real,
                           std::string ,
                           bool>())
        .add_property("filename", &HDF5File::getFilename,
                                  &HDF5File::setFilename)
		.add_property("iomode", &HDF5File::getIomode,
										&HDF5File::setIomode)
        .add_property("unfolded", &HDF5File::getUnfolded,
                                  &HDF5File::setUnfolded)
        .add_property("length_factor", &HDF5File::getLengthFactor,
                                       &HDF5File::setLengthFactor)
        .add_property("length_unit", &HDF5File::getLengthUnit,
                                     &HDF5File::setLengthUnit)
        .add_property("append", &HDF5File::getAppend,
                                  &HDF5File::setAppend)
        .def("write", &HDF5File::write)
      ;
    }
  }
}

