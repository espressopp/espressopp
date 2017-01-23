#  Copyright (c) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
*********************************
espressopp.analysis.SystemMonitor
*********************************

SystemMonitor prints and logs to file values obtained from Observables like
temperature, pressure or potential energy.

.. function:: espressopp.analysis.SystemMonitor(system, integrator, output)

            :param system: The system object.
            :type system: espressopp.System
            :param integrator: The MD integrator.
            :type integrator: espressopp.integrator.MDIntegrator
            :param output: The output object.
            :type output: espressopp.analysis.SystemMonitorOutputCSV

.. function:: espressopp.analysis.SystemMonitor.add_observable(name, observable, is_visible)

            The function adds new observable to SystemMonitor.

            :param name: The name of observable
            :type name: str
            :param observable: The observable, eg. espressopp.analysis.PotentialEnergy
            :type name:
            :param is_visible: If set to True then values will be print on console.
            :type is_visible: bool

.. function:: espressopp.analysis.SystemMonitor.info()

            The method print out on console the values of observables.

**CSV Output**

The output of SystemMonitor to CSV files.

.. function:: espressopp.analysis.SystemMonitorOutputCSV(file_name, delimiter)

            :param file_name: The name of CSV file.
            :type file_name: str
            :param delimiter: The field delimiter, by default it is tabulator.
            :type delimiter: str


Example

>>> interaction = espressopp.interaction.VerletListLennardJones(verletlist)
>>> interaction.setPotential(type1=0, type2=0,
                             potential=espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0,
                                                                           cutoff=2.0))
>>> system_monitor_csv = espressopp.analysis.SystemMonitorOutputCSV('out.csv')
>>> system_monitor = espressopp.analysis.SystemMonitor(
        system, integrator, espressopp.analysis.SystemMonitorOutputCSV('out.csv'))
>>> system_monitor.add_observable('pot', espressopp.analysis.PotentialEnergy(system, interaction))
>>> ext_analysis = espressopp.integrator.ExtAnalyze(system_monitor, 10)
>>> integrator.addExtension(ext_analysis)

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.AnalysisBase import *  #NOQA
from _espressopp import analysis_SystemMonitor
from _espressopp import analysis_SystemMonitorOutputCSV


class SystemMonitorOutputCSVLocal(analysis_SystemMonitorOutputCSV):
    def __init__(self, file_name, delimiter='\t'):
        if pmi.workerIsActive():
            cxxinit(self, analysis_SystemMonitorOutputCSV, file_name, delimiter)


class SystemMonitorLocal(analysis_SystemMonitor):
    def __init__(self, system, integrator, output):
        if pmi.workerIsActive():
            cxxinit(self, analysis_SystemMonitor, system, integrator, output)

    def add_observable(self, name, observable, is_visible=True):
        if pmi.workerIsActive():
            self.cxxclass.add_observable(self, name, observable, is_visible)

    def info(self):
        if pmi.workerIsActive():
            self.cxxclass.info(self)

    def dump(self):
        if pmi.workerIsActive():
            self.cxxclass.dump(self)

if pmi.isController:
    class SystemMonitorOutputCSV:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espressopp.analysis.SystemMonitorOutputCSVLocal')

    class SystemMonitor:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.SystemMonitorLocal',
            pmicall=['add_observable', 'info', 'dump']
            )
