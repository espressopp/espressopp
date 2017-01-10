/*
 Copyright (C) 2012-2016
     Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011
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

#ifndef _LBOUTPUT_SCREEN_HPP
#define _LBOUTPUT_SCREEN_HPP

#include "LBOutput.hpp"
#include "esutil/Timer.hpp"

namespace espressopp {
   namespace analysis {
      class LBOutputScreen : public LBOutput {
      public:
         LBOutputScreen(shared_ptr<System> _system,
                        shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann);

         void writeOutput();
         void findLBMom(int _mode);

         void setMDMom(Real3D _mdMom);
         Real3D getMDMom();

         void setLBMom(Real3D _lbMom);
         Real3D getLBMom();

         void setLBTimerOld(real _lbTimerOld);
         real getLBTimerOld();

         void setLBTimerNew(real _lbTimerNew);
         real getLBTimerNew();

         void setOldStepNum(long int _oldStepNum);
         long int getOldStepNum();

         static void registerPython();

      private:
         int oldStepNum;
         real lbTimerOld, lbTimerNew;
         Real3D lbMom, mdMom;
         esutil::WallTimer timeLBtoMD;  //!< used for timing
      };
   }
}

#endif
