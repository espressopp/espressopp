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

// ESPP_CLASS
#ifndef _ANALYSIS_MEANSQUAREDISPL_HPP
#define _ANALYSIS_MEANSQUAREDISPL_HPP

#include "ConfigsParticleDecomp.hpp"

namespace espressopp {
    namespace analysis {

        /*
         * Class derived from ConfigsParticleDecomp.
         * 
         * This implementation of mean square displacement calculation does not take into
         * account particle masses. It is correct if all the particles have equal masses only.
         * Otherwise it should be modified.
         */

        class MeanSquareDispl : public ConfigsParticleDecomp {
        public:

            MeanSquareDispl(shared_ptr<System> system) : ConfigsParticleDecomp(system) {
                // by default 
                setPrint_progress(true);
                key = "unfolded";
            }

            MeanSquareDispl(shared_ptr<System> system, int chainlength) :
                                ConfigsParticleDecomp(system, chainlength) {
                // by default 
                setPrint_progress(true);
                key = "unfolded";
            }

            ~MeanSquareDispl() {
            }

            virtual python::list compute() const;
            python::list computeG2() const;
            python::list computeG3() const;

            void setPrint_progress(bool _print_progress) {
                print_progress = _print_progress;
            }

            bool getPrint_progress() {
                return print_progress;
            }

            static void registerPython();
        private:
            bool print_progress;
            void printReal3D(Real3D v) const {
//                real x, y, z;
//                x = v.getItem(0); y = v.getItem(1); z = v.getItem(2);
//                printf("(%10.4f,%10.4f,%f10.4)", x, y, z);
                printf("(%10.4f,%10.4f,%f10.4)", v.getItem(0), v.getItem(1), v.getItem(2));
            };
        };
    }
}

#endif
