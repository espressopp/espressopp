/*
  Copyright (c) 2015-2016
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
#include "DumpTopology.hpp"
#include <vector>


namespace espressopp {
namespace io {

void DumpTopology::ObserveTuple(shared_ptr<FixedPairList> fpl) {
  fpls_.push_back(fpl);
}

void DumpTopology::perform_action() {
  Dump();
}

void DumpTopology::Dump() {
  // Format: <step1><pair_idx_0><size><pid1><pid2><pid3><pid4><pair_idx_1><size><pid...><step2>..
  int idx = 0;
  for (std::vector<shared_ptr<FixedPairList> >::iterator it = fpls_.begin();
       it != fpls_.end(); ++it) {
    fpl_buffer_.push_front(integrator_->getStep());
    std::vector<longint> pairs = (*it)->getPairList();
    fpl_buffer_.push_front(idx);
    fpl_buffer_.push_front(pairs.size()/2);
    for (std::vector<longint>::iterator itt = pairs.begin(); itt != pairs.end();) {
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
    }
    idx++;
  }
}

void DumpTopology::ClearBuffer() {
  fpl_buffer_.clear();
}

python::list DumpTopology::GetData() {
  python::list ret;
  for (FplBuffer::iterator it = fpl_buffer_.begin(); it != fpl_buffer_.end(); ++it) {
    ret.append(*it);
  }
  return ret;
}

// Python wrapping
void DumpTopology::registerPython() {
  using namespace espressopp::python;  //NOLINT

  class_<DumpTopology, bases<ParticleAccess>, boost::noncopyable>
      ("io_DumpTopology", init<shared_ptr<System>, shared_ptr<integrator::MDIntegrator> >())
      .def("clear_buffer", &DumpTopology::ClearBuffer)
      .def("get_data", &DumpTopology::GetData)
      .def("dump", &DumpTopology::Dump)
      .def("observe_tuple", &DumpTopology::ObserveTuple);
}

}  // namespace io
}  // namespace espressopp
