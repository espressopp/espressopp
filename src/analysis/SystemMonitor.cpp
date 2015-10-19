/*
  Copyright (c) 2015
    Jakub Krajniak (jkrajniak at gmail.com)

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

#include <string>
#include <utility>
#include <vector>
#include "python.hpp"
#include "SystemMonitor.hpp"
#include "integrator/MDIntegrator.hpp"

namespace espressopp{
namespace analysis {

void SystemMonitor::perform_action() {
    current_step_ = integrator_->getStep();
    values_.clear();
    values_.push_back(current_step_);
    values_.push_back(current_step_ * integrator_->getTimeStep());

    compute_kinetic_energy();
    compute_observables();

    write();
}

void SystemMonitor::compute_observables() {
  for (ObservableList::iterator it = observables_.begin(); it != observables_.end(); ++it) {
    values_.push_back(it->second->compute_real());
  }
}

void SystemMonitor::compute_kinetic_energy() {
  real T = temp_->computeRaw();
  real ekin = (3.0/2.0) * npart_->compute_real() * T;
  values_.push_back(T);
  values_.push_back(ekin);
}

void SystemMonitor::write() {
  // That's a text file, let's only write on one node.
  if (system_->comm->rank() == 0) {
    std::ofstream output_file;
    std::stringstream ss;
    // First run, write header.
    if (!header_written_) {  // First run;
      output_file.open(file_name_.c_str(), std::fstream::out);
      for (std::vector<std::string>::iterator it = header_.begin(); it != header_.end(); ++it) {
        ss << *it;
        if (it != header_.end()-1)
            ss << delimiter_;
      }
      ss << std::endl;;
      header_written_ = true;
    } else {
      output_file.open(file_name_.c_str(), std::ofstream::out | std::ofstream::app);
    }
    
    for (std::vector<real>::iterator it = values_.begin(); it != values_.end(); ++it) {
      ss << *it;
      if (it != values_.end()-1)
        ss << delimiter_;
    }
    ss << std::endl;

    // Close file.
    output_file << ss.str();
    output_file.close();
  }
}

void SystemMonitor::info() {
  if (system_->comm->rank() == 0) {
    int idx = 0;
    if (!header_shown_) {
      for (std::vector<std::string>::iterator it = header_.begin(); it != header_.end(); ++it) {
        if (visible_observables_[idx] == 1) {
          std::cout << *it;
          if (it != header_.end()-1)
              std::cout << "\t";
        }
        idx++;
      }
      std::cout << std::endl;
      header_shown_ = true;
    }
    // Print data
    idx = 0;
    for (std::vector<real>::iterator it = values_.begin(); it != values_.end(); ++it) {
      if (visible_observables_[idx] == 1) {
        std::cout << *it;
        if (it != values_.end()-1)
          std::cout << "\t";
      }
      idx++;
    }
    std::cout << std::endl;
  }
}

void SystemMonitor::add_observable(std::string name, shared_ptr<Observable> obs, 
    bool is_visible) {
  observables_.push_back(std::make_pair(name, obs));
  header_.push_back(name);
  if (is_visible)
    visible_observables_.push_back(1);
  else
    visible_observables_.push_back(0);
}

void SystemMonitor::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<SystemMonitor, bases< ParticleAccess > >
      ("analysis_SystemMonitor", init<
          shared_ptr<System>,
          shared_ptr<integrator::MDIntegrator>,
          std::string,
          std::string
          >())
      .def("add_observable", &SystemMonitor::add_observable)
      .def("info", &SystemMonitor::info)
      .def("dump", &SystemMonitor::perform_action);
}

/** Implementation of SystemMonitorOutputs. **/

SystemMonitorOutputCSV::write() {

}

}  // end namespace analysis
}  // end namespace espresso
