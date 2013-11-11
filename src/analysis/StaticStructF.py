from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_StaticStructF

class StaticStructFLocal(ObservableLocal, analysis_StaticStructF):
  'The (local) compute the static structure function.'
  def __init__(self, system):
    cxxinit(self, analysis_StaticStructF, system)
    
  def compute(self, nqx, nqy, nqz, bin_factor, ofile = None):
    if ofile is None:
      return self.cxxclass.compute(self, nqx, nqy, nqz, bin_factor)
    else:    
      #run compute on each CPU
      result = self.cxxclass.compute(self, nqx, nqy, nqz, bin_factor)
      #create the outfile only on CPU 0
      if pmi.isController:
        myofile = 'qsq_' + str(ofile) + '.txt'
        outfile = open (myofile, 'w')
        for i in range (len(result)):
          line = str(result[i][0]) + "\t" + str(result[i][1]) + "\n"
          outfile.write(line)
        outfile.close()
      return result   
   
if pmi.isController :
  class StaticStructF(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espresso.analysis.StaticStructFLocal'
    )
