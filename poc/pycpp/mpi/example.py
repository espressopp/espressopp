# LD_LIBRARY_PATH=$(BOOST_LIB)
import sys
sys.path.append('python')

import mpi_test
print sys.modules['mpi_test']

mpi_test.test(1)
