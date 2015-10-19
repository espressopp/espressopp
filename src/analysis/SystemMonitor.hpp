/*
  Copyright (c) 2015
      Jakub Krajniak (jkrajniak at gmail.com)
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

// ESPP_CLASS
#ifndef _ANALYSIS_SYSTEMENERGY_HPP
#define _ANALYSIS_SYSTEMENERGY_HPP

#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
#include <vector>
#include "types.hpp"
#include "integrator/MDIntegrator.hpp"
#include "ParticleAccess.hpp"
#include "Temperature.hpp"
#include "NPart.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "io/FileBackup.hpp"

namespace espressopp {
namespace analysis {


class SystemMonitorOutput {
 public:
  virtual void write() = 0;
}

class SystemMonitorOutputCSV: public SystemMonitorOutput {
 public:
  SystemMonitorOutputCSV(std::string file_name, std::string delimiter) :
      file_name_(file_name), delimiter_(delimiter) {
    header_written_ = false;
  }
  void write();

  static void registerPython();
 private:
  bool header_written_;
  std::string file_name_;
  std::string delimiter_;
};


class SystemMonitor : public ParticleAccess {
 public:
  typedef std::vector<std::pair<std::string, shared_ptr<Observable> > > ObservableList;
  SystemMonitor(shared_ptr< System > system,
                shared_ptr<integrator::MDIntegrator> integrator,
                shared_ptr<SystemMonitorOutput> output):
        ParticleAccess(system),
        system_(system),
        integrator_(integrator),
        output_(output) {
    temp_ = shared_ptr<Temperature>(new Temperature(system));
    npart_ = shared_ptr<NPart>(new NPart(system));

    header_shown_ = false;
    if (system->comm->rank() == 0) {
      header_.push_back("step");
      header_.push_back("time");
      header_.push_back("T");
      header_.push_back("Ekin");
      visible_observables_.push_back(1);
      visible_observables_.push_back(1);
      visible_observables_.push_back(1);
      visible_observables_.push_back(1);
    }
  }

  ~SystemMonitor() {
  }
  void perform_action();
  void info();

  static void registerPython();

 private:
  void write();
  void compute_observables();
  void compute_kinetic_energy();

  void add_observable(std::string name, shared_ptr<Observable> obs, bool is_visible);

  int current_step_;
  bool header_written_;
  bool header_shown_;
  std::vector<real> values_;
  std::vector<std::string> header_;
  std::vector<int> visible_observables_;
  shared_ptr<System> system_;
  shared_ptr<integrator::MDIntegrator> integrator_;
  shared_ptr<Temperature> temp_;
  shared_ptr<NPart> npart_;

  shared_ptr<SystemMonitorOutput> output_;
  ObservableList observables_;
};

}  // end namespace analysis
}  // end namespace espressopp
#endif
