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
#include "System.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "esutil/RNG.hpp"
#include "mpi.hpp"
#include "esutil/Error.hpp"

#include <limits>

#include <mpi4py/mpi4py.h>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
# define VT_ON()
# define VT_OFF()
#endif

namespace espressopp {

  System::System() {
    comm = mpiWorld;
    CommunicatorIsInitialized = false;
    
    maxCutoff = 0.0;
  }

  System::System(python::object _pyobj) {

    //This is defined in mpi4py.h
    //struct __pyx_obj_6mpi4py_3MPI_Comm {
    //  PyObject_HEAD
    //  MPI_Comm ob_mpi;
    //  int flags;
    //};
    
    // Following is some extreme typecasting which we need to convert the
    // pmi python object pmi._MPIcomm into a shared_ptr< boost::mpi::communicator >
    // I have not yet figured out how to do it in a more elegant way
    PyObject *pyobj = _pyobj.ptr();
    // in mpi4py.1.2.1 this has to be:
    // __pyx_obj_6mpi4py_3MPI_Comm * pyMPIComm = (__pyx_obj_6mpi4py_3MPI_Comm *) pyobj;
    // in version >= mpi4py.1.2.2 this has to be:
    PyMPICommObject* pyMPIComm = (PyMPICommObject*) pyobj;
    MPI_Comm * comm_p = &pyMPIComm->ob_mpi;
    shared_ptr< mpi::communicator > newcomm = make_shared< mpi::communicator >(*comm_p, mpi::comm_attach);

    comm = newcomm;
    maxCutoff = 0.0;
  }

  void System::setSkin(real _skin){
    skin = _skin;
    if(storage){
      // TODO !! It may give an error for decompositions different from 
      // DomainDecomposition
      // probably function getInt3DCellGrid() should be pure virtual in Storage
      Int3D cellGr = storage->getInt3DCellGrid();
      real cs = maxCutoff + skin;
      if( cs > std::min( std::min( cellGr[0], cellGr[1]), cellGr[2] ) ){
        storage -> cellAdjust();
      }
    }
    //storage -> decompose();  // it's not nessesary because at the end of cellAdjust()
                               // the signal onParticlesChanged is sent
  }
  real System::getSkin(){
    return skin;
  }

  void System::addInteraction(shared_ptr< interaction::Interaction > ia){
    shortRangeInteractions.push_back(ia);
    
    // check if the cutoff of this interaction is bigger then maxCutoff
    real cut = ia->getMaxCutoff();
    
    bool isCutoffInf = std::numeric_limits<real>::has_infinity &&
            cut == std::numeric_limits<real>::infinity();

    if( !isCutoffInf && cut > maxCutoff) maxCutoff=cut;
  }

  void System::removeInteraction(int i){
    esutil::Error err(comm);
    if(i>=shortRangeInteractions.size()){
      std::stringstream msg;
      msg << "Probably you are trying to remove the interaction "<<i<<
              " which does not exist. Check your script!";
      err.setException( msg.str() );
    }
    err.checkException();
    
    size_t iIter = i;
    real maxCutoffDelete = shortRangeInteractions[iIter]->getMaxCutoff();
    
    shortRangeInteractions.erase(shortRangeInteractions.begin()+i);
    
    // check if the maxCutoff changed or not
    if(maxCutoffDelete>=maxCutoff){
      maxCutoff=0.0;
      for (size_t j = 0; j < shortRangeInteractions.size(); j++) {
        real cut = shortRangeInteractions[j]->getMaxCutoff();
        maxCutoff = std::max(maxCutoff, cut);
      }
    }
  }

  shared_ptr< interaction::Interaction > System::getInteraction(int i)
  {
    return shortRangeInteractions[i];
  }

  int System::getNumberOfInteractions()
  {
	return shortRangeInteractions.size();
  }

  /* If one wants overload scaleVolume it should be done here as well:
   * - Storage
   * - BC
   * - DomainDecomposition, etc.
   * - CellGrid
   */
  
  // Scale all coordinates of the system, isotropic case (cubic box)
  void System::scaleVolume(real s, bool particleCoordinates) {
	// the size of the system should be modified first because of the cell size recalculation
    // in xDecomposition
    bc->scaleVolume(s);
	storage->scaleVolume(s, particleCoordinates);
  }

  // Scale all coordinates of the system, anisotropic case (rectangular system!!!).
  // Now the scale parameter is vector (Real3D). It is the case of the anisotropic extension of the
  // system. Certainly it should be modified for triclinic box. Scale parameter s should be not
  // Real3D but tensor.
  void System::scaleVolume(Real3D s, bool particleCoordinates) {
	// the size of the system should be modified first because of the cell size recalculation
    // in xDecomposition
	bc->scaleVolume(s);
	storage->scaleVolume(s, particleCoordinates);
  }
  
  void System::setTrace(bool flag) {
     if (flag){
        VT_ON();
     }
     else{
        VT_OFF();
     }
  }
  
  /////////////////////////////////////////////////////
  // Helper Function for Python interface  ////////////
  /////////////////////////////////////////////////////
  // this function is called from python by default
  void System::scaleVolume3D(Real3D s) {
    if(s[0]==s[1] && s[0]==s[2]){
      real s1 = s[0];
      scaleVolume(s1, true);
    }
    else
      scaleVolume(s, true);
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  System::registerPython() {
    using namespace espressopp::python;

    class_< System > ("System", init<>())
      .add_property("skin", &System::getSkin, &System::setSkin)
    
      .def(init< python::object >())
      .def_readwrite("storage", &System::storage)
      .def_readwrite("bc", &System::bc)
      .def_readwrite("rng", &System::rng)
//      .def_readwrite("shortRangeInteractions",
//		     &System::shortRangeInteractions)
      .def_readonly("maxCutoff", &System::maxCutoff)
      .def("addInteraction", &System::addInteraction)
      .def("removeInteraction", &System::removeInteraction)
      .def("getInteraction", &System::getInteraction)
      .def("getNumberOfInteractions", &System::getNumberOfInteractions)
      .def("scaleVolume", &System::scaleVolume3D)
      .def("setTrace", &System::setTrace)
      ;
  }
}
