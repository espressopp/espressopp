#  Copyright (C) 2016
#      Jakub Krajniak <jkrajniak at gmail.com>
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


r"""
*******************************
**espressopp.VerletListHybrid**
*******************************
"""
from espressopp import pmi
import _espressopp
from espressopp.esutil import cxxinit


class VerletListHybridATLocal(_espressopp.VerletListHybridAT):
    def __init__(self, system, cutoff, exclusionlist=None):

        if pmi.workerIsActive():
            if exclusionlist is None:
                # rebuild list in constructor
                cxxinit(self, _espressopp.VerletListHybridAT,
                        system, cutoff, True)
            else:
                # do not rebuild list in constructor
                cxxinit(self, _espressopp.VerletListHybridAT,
                        system, cutoff, False)
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
                # now rebuild list with exclusions
                self.cxxclass.rebuild(self)

    def totalSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)

    def localSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.localSize(self)

    def exclude(self, exclusionlist):
        """
        Each processor takes the broadcasted exclusion list
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pair in exclusionlist:
                pid1, pid2 = pair
                self.cxxclass.exclude(self, pid1, pid2)
            # rebuild list with exclusions
            self.cxxclass.rebuild(self)

    def getAllPairs(self):
        if pmi.workerIsActive():
            pairs = []
            npairs = self.localSize()
            for i in range(npairs):
                pair = self.cxxclass.getPair(self, i+1)
                pairs.append(pair)
            return pairs


class VerletListHybridCGLocal(_espressopp.VerletListHybridCG):
    def __init__(self, system, cutoff, exclusionlist=None):

        if pmi.workerIsActive():
            if exclusionlist is None:
                # rebuild list in constructor
                cxxinit(self, _espressopp.VerletListHybridCG,
                        system, cutoff, True)
            else:
                # do not rebuild list in constructor
                cxxinit(self, _espressopp.VerletListHybridCG,
                        system, cutoff, False)
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
                # now rebuild list with exclusions
                self.cxxclass.rebuild(self)

    def totalSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)

    def localSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.localSize(self)

    def exclude(self, exclusionlist):
        """
        Each processor takes the broadcasted exclusion list
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pair in exclusionlist:
                pid1, pid2 = pair
                self.cxxclass.exclude(self, pid1, pid2)
            # rebuild list with exclusions
            self.cxxclass.rebuild(self)

    def getAllPairs(self):
        if pmi.workerIsActive():
            pairs = []
            npairs = self.localSize()
            for i in range(npairs):
                pair = self.cxxclass.getPair(self, i+1)
                pairs.append(pair)
            return pairs


if pmi.isController:
    class VerletListHybridAT(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.VerletListHybridATLocal',
            pmiproperty=['builds'],
            pmicall=['totalSize', 'exclude', 'connect',
                     'disconnect', 'getVerletCutoff'],
            pmiinvoke=['getAllPairs']
        )

    class VerletListHybridCG(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.VerletListHybridCGLocal',
            pmiproperty=['builds'],
            pmicall=['totalSize', 'exclude', 'connect',
                     'disconnect', 'getVerletCutoff'],
            pmiinvoke=['getAllPairs']
        )
