# qwghlm
CPPFLAGS=\
	-I /home/olenz/software/include/boost-1_35 \
	-I /usr/include/mpi
LDLIBS=-L /home/olenz/software/lib/ \
	-lboost_mpi-gcc42-mt \
	-lboost_serialization-gcc42-mt \
	-lmpi
CXXFLAGS=-O3

# MPIP
#CPPFLAGS=\
#	-I /people/thnfs/homes/lenzo/software/include/boost-1_35 \
#	-I /people/thnfs/homes/lenzo/software/include/
#LDFLAGS=\
#	-L /people/thnfs/homes/lenzo/software/lib \
#	-lboost_mpi-gcc42-mt \
#	-lboost_serialization-gcc42-mt \
#	-lmpi
#CXXFLAGS=-O3


all: example

clean:
	-rm -f *.o example
