
from espressopp.integrator.SampleBase import *

class DerivedSampleBaseLocal(SampleBaseLocal):
    def __init__(self):
        if pmi.workerIsActive():
            print('DerivedSampleBaseLocal.__init__')

if pmi.isController:
    class DerivedSampleBase(SampleBase, metaclass=pmi.Proxy):
        pmiproxydefs = {
            'cls': 'espressopp.integrator.DerivedSampleBaseLocal'
        }
