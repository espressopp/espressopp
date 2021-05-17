/*
  Copyright (C) 2017,2020
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

#include "DumpXTCAdress.hpp"
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
#include <string>

#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"
#include "analysis/ConfigurationsExtAdress.hpp"

namespace espressopp
{
namespace io
{
bool DumpXTCAdress::open(const char *mode)
{
    if (mode[0] == 'a' && !boost::filesystem::exists(file_name))
    {
        // open file in append mode "w"
        fio = open_trx(file_name.c_str(), "w");
    }
    else
    {
        fio = open_trx(file_name.c_str(), mode);
    }

    return true;
}

void DumpXTCAdress::close() { close_trx(fio); }

void DumpXTCAdress::dump()
{
    std::shared_ptr<System> system = getSystem();
    analysis::ConfigurationsExtAdress conf(system, ftpl);
    conf.setUnfolded(unfolded);
    conf.gather();

    if (system->comm->rank() == 0)
    {
        t_trxframe frame;

        analysis::ConfigurationExtPtr conf_real = conf.back();

        int num_of_particles = conf_real->getSize();

        if (this->open("a"))
        {
            analysis::ConfigurationExtIterator cei = conf_real->getIterator();
            matrix box;
            rvec *coord = new rvec[num_of_particles];
            RealND props;
            props.setDimension(cei.currentProperties().getDimension());

            for (int i = 0; i < num_of_particles; i++)
            {
                props = cei.nextProperties();
                coord[i][0] = props[0] * length_factor;  // We only write coordinates to .xtc
                coord[i][1] = props[1] * length_factor;
                coord[i][2] = props[2] * length_factor;
            }

            int step = integrator->getStep();
            float time = integrator->getTimeStep() * step;

            frame.natoms = num_of_particles;
            frame.bTime = true;
            frame.time = time;
            frame.bStep = true;
            frame.step = step;
            frame.x = coord;
            frame.bLambda = false;
            frame.bAtoms = false;
            frame.bPrec = true;
            frame.prec = xtcprec;
            frame.bX = true;
            frame.bF = false;
            frame.bBox = true;

            // TODO(anyone): this works only for orthorombic BC
            Real3D bl = system->bc->getBoxL();

            for (int i = 0; i < dim; i++)
            {
                frame.box[i][0] = 0.;
                frame.box[i][1] = 0.;
                frame.box[i][2] = 0.;

                frame.box[i][i] = bl[i];
            }

            write_trxframe(fio, &frame, nullptr);

            delete[] coord;

            this->close();
        }
        else
        {
            std::cout << "Unable to open file: " << file_name << std::endl;
        }
    }
}

// Python wrapping
void DumpXTCAdress::registerPython()
{
    using namespace espressopp::python;  // NOLINT

    class_<DumpXTCAdress, bases<ParticleAccess>, boost::noncopyable>(
        "io_DumpXTCAdress",
        init<std::shared_ptr<System>, std::shared_ptr<FixedTupleListAdress>,
             std::shared_ptr<integrator::MDIntegrator>, std::string, bool, real, bool>())
        .add_property("filename", &DumpXTCAdress::getFilename, &DumpXTCAdress::setFilename)
        .add_property("unfolded", &DumpXTCAdress::getUnfolded, &DumpXTCAdress::setUnfolded)
        .add_property("append", &DumpXTCAdress::getAppend, &DumpXTCAdress::setAppend)
        .def("dump", &DumpXTCAdress::dump);
}
}  // end namespace io
}  // end namespace espressopp
