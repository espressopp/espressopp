/*
  Copyright (c) 2015-2017
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
#include <boost/mpi/timer.hpp>

namespace espressopp {
namespace analysis {

/**
 * Base class for SystemMonitor output formats.
 */
class SystemMonitorOutput {
  friend class SystemMonitor;
  public:
    virtual void write() = 0;
    void setSystem(shared_ptr<System> system) {
      system_ = system;
    }

    void setKeys(shared_ptr<std::vector<std::string> > keys) {
      keys_ = keys;
    }

    void setValues(shared_ptr<std::vector<real> > val) {
      values_ = val;
    }

    static void registerPython();

  protected:
    shared_ptr<std::vector<std::string> > keys_;
    shared_ptr<std::vector<real> > values_;
    shared_ptr<System> system_;
    static LOG4ESPP_DECL_LOGGER(theLogger);
};

/**
 * Output of SystemMonitor writen to CSV file.
 *
 * @param file_name The output file name.
 * @param delimiter The separator of fields.
 */
class SystemMonitorOutputCSV : public SystemMonitorOutput {
 public:
  SystemMonitorOutputCSV(std::string file_name, std::string delimiter) :
      file_name_(file_name), delimiter_(delimiter) {
    header_written_ = false;
  }
  void write();

  static void registerPython();

 private:
  std::string file_name_;
  std::string delimiter_;
  bool header_written_;
};

class SystemMonitorOutputDummy : public SystemMonitorOutput {
 public:
  void write() { }
  static void registerPython();
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
    header_ = make_shared<std::vector<std::string> >();
    values_ = make_shared<std::vector<real> >();

    output_->setKeys(header_);
    output_->setValues(values_);
    output_->setSystem(system);

    header_shown_ = false;
    if (system->comm->rank() == 0) {
      header_->push_back("step");
      header_->push_back("time");
      visible_observables_.push_back(1);
      visible_observables_.push_back(1);
    }
    elapsed_time_ = true;
  }

  ~SystemMonitor() { }

  void perform_action();
  void info();

  static void registerPython();

 private:
  void computeObservables();

  void addObservable(std::string name, shared_ptr<Observable> obs, bool is_visible);

  int current_step_;
  bool header_shown_;
  bool elapsed_time_;
  boost::mpi::timer timer_;

  shared_ptr<std::vector<real> > values_;
  shared_ptr<std::vector<std::string> > header_;
  std::vector<int> visible_observables_;
  shared_ptr<System> system_;
  shared_ptr<integrator::MDIntegrator> integrator_;

  shared_ptr<SystemMonitorOutput> output_;
  ObservableList observables_;

  real total_energy_;
  real potential_energy_;
};

}  // end namespace analysis
}  // end namespace espressopp
#endif
