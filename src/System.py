from espresso import pmi
from espresso.esutil import cxxinit

import _espresso
import MPI

class SystemLocal(_espresso.System):
    'The (local) System.'
    def __init__(self, commoid, pmicomm=None):
        'Local construction of a System'
        lcomm = None
        if pmicomm is None :
            lcomm=pmi._backtranslateOID(commoid)
        else :
            lcomm=pmicomm
#        print "CPU%d: System initialized on" % pmi._MPIcomm.rank,lcomm.getMPIcpugroup()
        comm = lcomm.getMPIsubcomm
        if comm != MPI.COMM_NULL :
            print "CPU%d: System initialized" % pmi._MPIcomm.rank
            cxxinit(self, _espresso.System, comm)

    def addInteraction(self, interaction):
        'add a short range list interaction'
        return self.cxxclass.addInteraction(self, interaction)

def getpmioid(obj):
    return obj.__pmioid

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin', 'shortRangeInteractions' ],
            pmicall = ['addInteraction' ]
            )

        def __init__(self, pmicomm) :
            self.pmiobjectclassdef=SystemLocal

        def __call__(self, method_self, *args, **kwds):
            method_self.pmiobjectclassdef = self.pmiobjectclassdef
            pmiobjectclass = pmi._translateClass(self.pmiobjectclassdef)
            pmicommOID = getpmioid(args[0].localcomm)
            method_self.pmiobject = pmi.create(pmiobjectclass, pmicommOID, __pmictr_pmicomm=pmicomm.localcomm)
            method_self.pmiobject._pmiproxy = method_self

#        def __init__(self, pmiobjectclassdef):
#            self.pmiobjectclassdef = pmiobjectclassdef
#        def __call__(self, method_self, *args, **kwds):
#            method_self.pmiobjectclassdef = self.pmiobjectclassdef
#            pmiobjectclass = _translateClass(self.pmiobjectclassdef)
#            method_self.pmiobject = create(pmiobjectclass, *args, **kwds)
#            method_self.pmiobject._pmiproxy = method_self
