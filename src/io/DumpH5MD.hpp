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

// ESPP_CLASS
#ifndef _IO_DUMPH5MD_HPP
#define _IO_DUMPH5MD_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "FixedTupleListAdress.hpp"
#include <Python.h>
#include <object.h>

namespace espressopp {
namespace io {

template <class T> void init_pb(Py_buffer *pb, int ndim, int *shape);
void free_pb(Py_buffer *pb);

/** Store particle data in Python buffers */
class DumpH5MD : public SystemAccess {
 public:
  DumpH5MD(shared_ptr<System> system, bool is_adress);
  ~DumpH5MD();

  int get_NLocal() { return NLocal; }
  void set_store_position(bool _s) { store_position = _s; }
  bool get_store_position() { return store_position; }
  void set_store_species(bool _s) { store_species = _s; }
  bool get_store_id() { return store_id; }
  void set_store_id(bool _s) { store_id = _s; }
  bool get_store_species() { return store_species; }
  void set_store_state(bool _s) { store_state = _s; }
  bool get_store_state() { return store_state; }
  void set_store_velocity(bool _s) { store_velocity = _s;}
  bool get_store_velocity() { return store_velocity; }
  void set_store_force(bool _s) { store_force = _s; }
  bool get_store_force() { return store_force; }
  void set_store_charge(bool _s) { store_charge = _s; }
  bool get_store_charge() { return store_charge; }
  void set_store_lambda(bool _s) { store_lambda = _s;}
  bool get_store_lambda() { return store_lambda; }
  void set_store_res_id(bool _s) { store_res_id = _s; }
  bool get_store_res_id() { return store_res_id; }
  void set_store_mass(bool _s) { store_mass = _s; }
  bool get_store_mass() { return store_mass; }

  void clear_buffers();

  void update();
  PyObject* getPosition();
  PyObject* getId();
  PyObject* getImage();
  PyObject* getSpecies();
  PyObject* getState();
  PyObject* getVelocity();
  PyObject* getForce();
  PyObject* getMass();
  PyObject* getCharge();
  PyObject* getLambda();
  PyObject* getResId();

  static void registerPython();

 private:
  bool store_position, store_velocity, store_mass, store_id, store_force;
  bool store_species, store_state, store_charge, store_lambda;
  bool store_res_id;

  bool is_adress_;
  bool cleared;
  int NLocal;
  Py_buffer position;
  Py_buffer image;
  Py_buffer velocity;
  Py_buffer mass;
  Py_buffer id;
  Py_buffer force;
  Py_buffer species;
  Py_buffer state;
  Py_buffer charge;
  Py_buffer lambda;
  Py_buffer res_id;
};
}  // namespace io
}  // namespace espressopp

#endif
