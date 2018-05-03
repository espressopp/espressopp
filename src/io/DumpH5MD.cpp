/*
  Copyright (c) 2015-2018
      Jakub Krajniak (jkrajniak at gmail.com)

  Copyright (c) 2015
      Pierre de Buyl

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

#include <Python.h>
#include <iostream>
#include <vector>
#include "DumpH5MD.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

const char *get_format(double x) { return "d"; }
const char *get_format(float x) { return "f"; }
const char *get_format(int x) { return "i"; }

namespace espressopp {
  namespace io {

    /** Initialize a Py_buffer object
     */
    template <class T> void init_pb(Py_buffer *pb, int ndim, int *shape)
    {
      T dum=0;
      // Setting basic variable of the Py_buffer
      pb->suboffsets = NULL;
      pb->internal = NULL;
      pb->obj = NULL;
      pb->readonly = 1;
      pb->ndim = ndim;
      // The format is computed as a function of the template type T
      pb->format = (char*)malloc(2*sizeof(const char));
      strcpy(pb->format, get_format(dum));
      // Allocation and setting of the shape, stride and len
      pb->shape = (Py_ssize_t *)malloc(pb->ndim*sizeof(Py_ssize_t));
      pb->strides = NULL;
      pb->itemsize = sizeof(T);
      int old_len = pb->len;
      pb->len = pb->itemsize;
      for (int i=0; i<pb->ndim; i++) {
        pb->shape[i] = shape[i];
        pb->len *= pb->shape[i];
      }
      pb->buf = malloc(pb->len);
    }

    void free_pb(Py_buffer *pb) {
      free(pb->format);
      pb->format = NULL;
      free(pb->shape);
      pb->shape = NULL;
      free(pb->buf);
      pb->buf = NULL;
      pb->len = 0;
    }

    DumpH5MD::DumpH5MD(shared_ptr<System> system, bool is_adress) :
        SystemAccess(system), is_adress_(is_adress) {

      store_position = true;
      store_velocity = store_mass = store_force = store_species = store_state = false;
      store_lambda = false;
      store_res_id = false;

      cleared = true;
      NLocal = -1;
      position.len = 0;
      velocity.len = 0;
      mass.len = 0;
      id.len = 0;
      force.len = 0;
      species.len = 0;
      state.len = 0;
      lambda.len = 0;
      res_id.len = 0;
    }

    DumpH5MD::~DumpH5MD() {
      clear_buffers();
    }

    void DumpH5MD::clear_buffers() {
      if (!cleared) {
        if (store_position) {
          free_pb(&position);
          free_pb(&image);
        }
        free_pb(&id);
        if (store_species) free_pb(&species);
        if (store_state) free_pb(&state);
        if (store_velocity) free_pb(&velocity);
        if (store_charge) free_pb(&charge);
        if (store_force) free_pb(&force);
        if (store_lambda) free_pb(&lambda);
        if (store_res_id) free_pb(&res_id);
        if (store_mass) free_pb(&mass);
        cleared = true;
      }
    }

    void DumpH5MD::update() {
      System& system = getSystemRef();
      int shape[2];

      clear_buffers();

      NLocal = system.storage->getNRealParticles();

      shape[0] = NLocal;
      shape[1] = 3;

      if (store_position) {
        init_pb<real>(&position, 2, shape);
        init_pb<longint>(&image, 2, shape);
      }
      init_pb<longint>(&id, 1, shape);
      if (store_mass) init_pb<real>(&mass, 1, shape);
      if (store_species) init_pb<int>(&species, 1, shape);
      if (store_state) init_pb<int>(&state, 1, shape);
      if (store_velocity) init_pb<real>(&velocity, 2, shape);
      if (store_force) init_pb<real>(&force, 2, shape);
      if (store_charge) init_pb<real>(&charge, 1, shape);
      if (store_lambda) init_pb<real>(&lambda, 1, shape);
      if (store_res_id) init_pb<int>(&res_id, 1, shape);
      cleared = false;

      CellList realCells = system.storage->getRealCells();

      int i = 0;
      if (is_adress_) {
        shared_ptr<FixedTupleListAdress> ftpl = system.storage->getFixedTuples();
        FixedTupleListAdress::iterator it3;
        for (iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          it3 = ftpl->find(&(*cit));
          if (it3 != ftpl->end()) {
            std::vector<Particle*> atList = it3->second;
            for (std::vector<Particle*>::iterator itv = atList.begin();
                 itv != atList.end(); ++itv) {
              Particle &at = **itv;
              if (store_position) {
                Real3D &tmpPos = at.position();
                Int3D &tmpImage = at.image();
                ((real *) position.buf)[3*i] = tmpPos[0];
                ((real *) position.buf)[3*i+1] = tmpPos[1];
                ((real *) position.buf)[3*i+2] = tmpPos[2];
                // Store image
                ((longint *) image.buf)[3*i] = tmpImage[0];
                ((longint *) image.buf)[3*i+1] = tmpImage[1];
                ((longint *) image.buf)[3*i+2] = tmpImage[2];
              }

              ((longint *) id.buf)[i] = at.id();
              if (store_mass) ((real *) mass.buf)[i] = at.mass();
              if (store_species) ((int *) species.buf)[i] = at.type();
              if (store_state) ((int *) state.buf)[i] = at.state();
              if (store_velocity) {
                Real3D &tmpVel = at.velocity();
                ((real *) velocity.buf)[3*i] = tmpVel[0];
                ((real *) velocity.buf)[3*i+1] = tmpVel[1];
                ((real *) velocity.buf)[3*i+2] = tmpVel[2];
              }
              if (store_force) {
                Real3D &tmpF = at.force();
                ((real *) force.buf)[3*i] = tmpF[0];
                ((real *) force.buf)[3*i+1] = tmpF[1];
                ((real *) force.buf)[3*i+2] = tmpF[2];
              }
              if (store_charge) ((real*) charge.buf)[i] = at.q();
              if (store_lambda) ((real*) lambda.buf)[i] = at.lambda();
              if (store_res_id) ((int*) res_id.buf)[i] = at.res_id();
              i++;
            }
          }
        }
      } else {
        for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          if (store_position) {
            Real3D &tmpPos = cit->position();
            Int3D &tmpImage = cit->image();
            ((real *) position.buf)[3*i] = tmpPos[0];
            ((real *) position.buf)[3*i+1] = tmpPos[1];
            ((real *) position.buf)[3*i+2] = tmpPos[2];
            // Store image
            ((longint *) image.buf)[3*i] = tmpImage[0];
            ((longint *) image.buf)[3*i+1] = tmpImage[1];
            ((longint *) image.buf)[3*i+2] = tmpImage[2];
          }
          ((longint *) id.buf)[i] = cit->getId();
          if (store_mass) ((real *) mass.buf)[i] = cit->getMass();
          if (store_species) ((int *) species.buf)[i] = cit->getType();
          if (store_state) ((int *) state.buf)[i] = cit->getState();
          if (store_velocity) {
            Real3D &tmpVel = cit->velocity();
            ((real *) velocity.buf)[3*i] = tmpVel[0];
            ((real *) velocity.buf)[3*i+1] = tmpVel[1];
            ((real *) velocity.buf)[3*i+2] = tmpVel[2];
          }
          if (store_force) {
            Real3D &tmpF = cit->force();
            ((real *) force.buf)[3*i] = tmpF[0];
            ((real *) force.buf)[3*i+1] = tmpF[1];
            ((real *) force.buf)[3*i+2] = tmpF[2];
          }
          if (store_charge) ((real*) charge.buf)[i] = cit->q();
          if (store_lambda) ((real*) lambda.buf)[i] = cit->lambda();
          if (store_res_id) ((int*) res_id.buf)[i] = cit->res_id();
          i++;
        }
      }
    }

    PyObject* DumpH5MD::getPosition() {
      if (store_position && position.len) 
        return PyMemoryView_FromBuffer(&position);
      Py_INCREF(Py_None);
      return Py_None;
    }
    PyObject* DumpH5MD::getImage() {
      if (store_position && image.len)
        return PyMemoryView_FromBuffer(&image);
      std::cout << "store_position=" << store_position << " position.len=" << position.len << " image.len=" << image.len << std::endl;
      std::cout << "WTF?" << std::endl;
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getId() {
      if (id.len) return PyMemoryView_FromBuffer(&id);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getSpecies() {
      if (store_species && species.len) return PyMemoryView_FromBuffer(&species);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getState() {
      if (store_state && state.len) return PyMemoryView_FromBuffer(&state);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getVelocity() {
      if (store_velocity && velocity.len)
        return PyMemoryView_FromBuffer(&velocity);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getForce() {
      if (store_force && force.len)
        return PyMemoryView_FromBuffer(&force);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getMass() {
      if (store_mass && mass.len)
        return PyMemoryView_FromBuffer(&mass);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getCharge() {
      if (store_charge && charge.len)
        return PyMemoryView_FromBuffer(&charge);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getLambda() {
      if (store_lambda && lambda.len)
        return PyMemoryView_FromBuffer(&lambda);
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* DumpH5MD::getResId() {
      if (store_res_id && res_id.len)
        return PyMemoryView_FromBuffer(&res_id);
      Py_INCREF(Py_None);
      return Py_None;
    }

    void DumpH5MD::registerPython() {
      using namespace espressopp::python;

      class_<DumpH5MD>
        ("io_DumpH5MD", init< shared_ptr< System >, bool >())
        .def("update", &DumpH5MD::update)
        .def("clear_buffers", &DumpH5MD::clear_buffers)
        .def("getPosition", &DumpH5MD::getPosition)
        .def("getImage", &DumpH5MD::getImage)
        .def("getId", &DumpH5MD::getId)
        .def("getSpecies", &DumpH5MD::getSpecies)
        .def("getState", &DumpH5MD::getState)
        .def("getVelocity", &DumpH5MD::getVelocity)
        .def("getForce", &DumpH5MD::getForce)
        .def("getMass", &DumpH5MD::getMass)
        .def("getCharge", &DumpH5MD::getCharge)
        .def("getLambda", &DumpH5MD::getLambda)
        .def("getResId", &DumpH5MD::getResId)
        .add_property("NLocal", &DumpH5MD::get_NLocal)
        .add_property("store_position", &DumpH5MD::get_store_position, &DumpH5MD::set_store_position)
        .add_property("store_species", &DumpH5MD::get_store_species, &DumpH5MD::set_store_species)
        .add_property("store_state", &DumpH5MD::get_store_state, &DumpH5MD::set_store_state)
        .add_property("store_velocity", &DumpH5MD::get_store_velocity, &DumpH5MD::set_store_velocity)
        .add_property("store_force", &DumpH5MD::get_store_force, &DumpH5MD::set_store_force)
        .add_property("store_charge", &DumpH5MD::get_store_charge, &DumpH5MD::set_store_charge)
        .add_property("store_lambda", &DumpH5MD::get_store_lambda, &DumpH5MD::set_store_lambda)
        .add_property("store_res_id", &DumpH5MD::get_store_res_id, &DumpH5MD::set_store_res_id)
        .add_property("store_mass", &DumpH5MD::get_store_mass, &DumpH5MD::set_store_mass)
        ;
    }
  }
}
