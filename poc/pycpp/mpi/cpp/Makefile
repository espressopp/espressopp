# To be able to run the examples, put 
# 

# MPIP:
# where to find the include files for boost and python
#CPPFLAGS=\
#	-I /people/thnfs/homes/lenzo/software/include/python2.5 \
#	-I /people/thnfs/homes/lenzo/software/include/boost-1_35
#
## where to find the boost and python libraries, plus their names
#LDFLAGS=-L /people/thnfs/homes/lenzo/software/lib \
#	-lboost_python-gcc42-mt -lpython2.5

# qwghlm:
CPPFLAGS=\
	-I /home/axel/boost/include/boost-1_35\
	-I /usr/include/python2.5

LIBS=\
	-L /home/axel/boost/lib \
	-lboost_mpi-gcc42-mt-1_35 \
	-lboost_serialization-gcc42-mt-1_35

PYLIBS=\
	-lboost_python-gcc42-mt-1_35 \
	-lpython2.5

# there is no need to change anything below this line
CXX=g++
MPXX=mpiCC
CC=gcc
CXXFLAGS=-fPIC -g -Wall
CFLAGS=-fPIC -g -Wall

all: mpi_test mpi_test.so mpi_test_embedded
clean:
	-rm -f *.o *.so
	-rm -f mpi_test mpi_test_embedded

mpi_test: main.o mpi_test_internal.o
	$(MPXX) -o $@ $^ $(LIBS)

mpi_test.so: mpi_test.o mpi_test_internal.o
	$(CXX) -shared $(CXXFLAGS) $(PYCPPFLAGS) $(LIBS) $(PYLIBS) -o $@ $^

mpi_test_embedded: mpi_test.o mpi_test_internal.o
	$(CXX) $(CXXFLAGS) $(PYCPPFLAGS) $(LIBS) $(PYLIBS) -o $@ $^
