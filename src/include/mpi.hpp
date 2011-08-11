#ifndef _MPI_HPP
#define _MPI_HPP

#include <boost/mpi.hpp>
#include <boost/smart_ptr.hpp>

extern boost::shared_ptr< boost::mpi::communicator > mpiWorld;

#endif
