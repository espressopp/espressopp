/*
  Copyright (c) 2015-2017,2021
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

#include "python.hpp"
#include <string>
#include <utility>
#include <vector>
#include "SystemMonitor.hpp"
#include "integrator/MDIntegrator.hpp"
#include <boost/format.hpp>

namespace espressopp
{
namespace analysis
{
void SystemMonitor::perform_action()
{
    current_step_ = integrator_->getStep();
    values_->clear();
    values_->push_back(current_step_);
    values_->push_back(current_step_ * integrator_->getTimeStep());

    computeObservables();
    if (system_->comm->rank() == 0)
    {
        output_->write();
    }
}

void SystemMonitor::computeObservables()
{
    total_energy_ = 0.0;
    potential_energy_ = 0.0;
    for (ObservableList::iterator it = observables_.begin(); it != observables_.end(); ++it)
    {
        int result_type = it->second->getResultType();
        if (result_type == Observable::real_vector)
        {
            std::vector<real> obs = it->second->compute_real_vector();
            for (int n = 0; n < it->second->getResultVectorSize(); n++)
            {
                values_->push_back(obs[n]);
            }
        }
        else if (result_type == Observable::real_scalar)
        {
            real val = it->second->compute_real();
            Observable::ObservableTypes obs_type = it->second->getObservableType();
            if (obs_type == Observable::POTENTIAL_ENERGY)
            {
                potential_energy_ += val;
            }
            else if (obs_type == Observable::KINETIC_ENERGY)
            {
                total_energy_ += val;
            }
            values_->push_back(val);
        }
    }
    total_energy_ += potential_energy_;
}

void SystemMonitor::info()
{
    if (system_->comm->rank() == 0)
    {
        int idx = 0;
        if (!header_shown_)
        {
            if (elapsed_time_)
            {
                std::cout << "elapsed\t";
            }
            for (std::vector<std::string>::iterator it = header_->begin(); it != header_->end();
                 ++it)
            {
                if (visible_observables_[idx] == 1)
                {
                    std::cout << *it;
                    if (it != header_->end() - 1) std::cout << "\t";
                }
                idx++;
            }
            std::cout << std::endl;
            header_shown_ = true;
        }
        // Print data
        idx = 0;
        if (elapsed_time_)
        {
            std::cout << timer_.elapsed() << "\t";
        }
        for (std::vector<real>::iterator it = values_->begin(); it != values_->end(); ++it)
        {
            if (visible_observables_[idx] == 1)
            {
                std::cout << *it;
                if (it != values_->end() - 1) std::cout << "\t";
            }
            idx++;
        }
        std::cout << std::endl;
    }
}

void SystemMonitor::addObservable(std::string name, shared_ptr<Observable> obs, bool is_visible)
{
    observables_.push_back(std::make_pair(name, obs));
    if (obs->getResultType() == Observable::real_scalar ||
        obs->getResultType() == Observable::old_format)
    {
        header_->push_back(name);
        if (is_visible)
        {
            visible_observables_.push_back(1);
        }
        else
        {
            visible_observables_.push_back(0);
        }
    }
    else if (obs->getResultType() == Observable::real_vector ||
             obs->getResultType() == Observable::int_vector)
    {
        boost::format frmt("%s[%d]");
        for (longint n = 0; n < obs->getResultVectorSize(); n++)
        {
            frmt.clear();
            frmt % name % n;
            header_->push_back(frmt.str());
            if (is_visible)
            {
                visible_observables_.push_back(1);
            }
            else
            {
                visible_observables_.push_back(0);
            }
        }
    }
}

void SystemMonitor::registerPython()
{
    using namespace espressopp::python;  // NOLINT
    class_<SystemMonitor, bases<ParticleAccess> >(
        "analysis_SystemMonitor", init<shared_ptr<System>, shared_ptr<integrator::MDIntegrator>,
                                       shared_ptr<SystemMonitorOutputCSV> >())
        .add_property("total_energy", make_getter(&SystemMonitor::total_energy_))
        .add_property("potential_energy", make_getter(&SystemMonitor::potential_energy_))
        .def("add_observable", &SystemMonitor::addObservable)
        .def("info", &SystemMonitor::info)
        .def("dump", &SystemMonitor::perform_action);
}

/** Implementation of SystemMonitorOutputs. **/

void SystemMonitorOutput::registerPython()
{
    using namespace espressopp::python;  // NOLINT
    class_<SystemMonitorOutput, boost::noncopyable>("analysis_SystemMonitorOutput", no_init);
}

void SystemMonitorOutputCSV::registerPython()
{
    using namespace espressopp::python;  // NOLINT
    class_<SystemMonitorOutputCSV, bases<SystemMonitorOutput> >("analysis_SystemMonitorOutputCSV",
                                                                init<std::string,  // file_name
                                                                     std::string   // delimiter
                                                                     >());
}

void SystemMonitorOutputDummy::registerPython()
{
    using namespace espressopp::python;  // NOLINT
    class_<SystemMonitorOutputDummy, bases<SystemMonitorOutput> >(
        "analysis_SystemMonitorOutputDummy", init<>());
}

void SystemMonitorOutputCSV::write()
{
    if (system_->comm->rank() == 0)
    {
        std::ofstream output_file;
        std::stringstream ss;
        if (!header_written_)
        {  // First run, write header;
            output_file.open(file_name_.c_str(), std::fstream::out);
            for (std::vector<std::string>::iterator it = keys_->begin(); it != keys_->end(); ++it)
            {
                ss << *it;
                if (it != keys_->end() - 1) ss << delimiter_;
            }
            ss << std::endl;
            ;
            header_written_ = true;
        }
        else
        {
            output_file.open(file_name_.c_str(), std::ofstream::out | std::ofstream::app);
        }

        // Write values;
        for (std::vector<real>::iterator it = values_->begin(); it != values_->end(); ++it)
        {
            ss << *it;
            if (it != values_->end() - 1) ss << delimiter_;
        }
        ss << std::endl;

        // Close file.
        output_file << ss.str();
        output_file.close();
    }
}
}  // end namespace analysis
}  // end namespace espressopp
