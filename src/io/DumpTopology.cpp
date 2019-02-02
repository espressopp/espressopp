/*
  Copyright (c) 2015-2018
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

void DumpTopology::observeTuple(shared_ptr<FixedPairList> fpl) {
  fpls_.push_back(fpl);
}

void DumpTopology::observeTriple(shared_ptr<FixedTripleList> ftl) {
  ftls_.push_back(ftl);
}

void DumpTopology::observeQuadruple(shared_ptr<FixedQuadrupleList> fql) {
  fqls_.push_back(fql);
}

void DumpTopology::saveDataToBuffer() {
  // Format: <step1><pair_idx_0><size><pid1><pid2><pid3><pid4><pair_idx_1><size><pid...><step2>..
  int idx = 0;
  for (std::vector<shared_ptr<FixedPairList> >::iterator it = fpls_.begin(); it != fpls_.end(); ++it) {
    std::vector<longint> pairs = (*it)->getPairList();
    fpl_buffer_.push_front(integrator_->getStep());
    fpl_buffer_.push_front(idx);
    fpl_buffer_.push_front(2);
    fpl_buffer_.push_front(pairs.size()/2);
    for (std::vector<longint>::iterator itt = pairs.begin(); itt != pairs.end();) {
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
    }
    idx++;
  }

  // Save data from triplet.
  idx = 0;
  for (std::vector<shared_ptr<FixedTripleList> >::iterator it = ftls_.begin(); it != ftls_.end(); ++it) {
    std::vector<longint> triples = (*it)->getTripleList();
    fpl_buffer_.push_front(integrator_->getStep());
    fpl_buffer_.push_front(idx);
    fpl_buffer_.push_front(3);
    fpl_buffer_.push_front(triples.size()/3);
    for (std::vector<longint>::iterator itt = triples.begin(); itt != triples.end();) {
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
    }
    idx++;
  }

  // Save data from quadruplets
  idx = 0;
  for (std::vector<shared_ptr<FixedQuadrupleList> >::iterator it = fqls_.begin(); it != fqls_.end(); ++it) {
    std::vector<longint> quadruples = (*it)->getQuadrupleList();
    fpl_buffer_.push_front(integrator_->getStep());
    fpl_buffer_.push_front(idx);
    fpl_buffer_.push_front(4);
    fpl_buffer_.push_front(quadruples.size()/4);
    for (std::vector<longint>::iterator itt = quadruples.begin(); itt != quadruples.end();) {
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
      fpl_buffer_.push_front(*(itt++));
    }
    idx++;
  }
}

void DumpTopology::clearBuffer() {
  fpl_buffer_.clear();
}

python::list DumpTopology::getData() {
  // TODO(jakub): direct access to vector from Python instead of rewriting it to boost Python list.
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
      .def("clear_buffer", &DumpTopology::clearBuffer)
      .def("get_data", &DumpTopology::getData)
      .def("dump", &DumpTopology::saveDataToBuffer)
      .def("observe_tuple", &DumpTopology::observeTuple)
      .def("observe_triple", &DumpTopology::observeTriple)
      .def("observe_quadruple", &DumpTopology::observeQuadruple);
}

}  // namespace io
}  // namespace espressopp
