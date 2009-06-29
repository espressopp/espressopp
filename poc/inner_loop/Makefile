CXXFLAGS=-O3

name=`hostname`

.PHONY: all clean benchmark

all: pairlists virtual_function iterators puetzloop

clean:
	-rm -f pairlists virtual_function iterators

%.s: %.C
	$(CXX) $(CXXFLAGS) -S $< -o $@

iterators: iterators_virtual.H iterators_nonvirtual.H iterators_block.H iterators_fobject.H

benchmark:
	./iterators | tee  $(name).benchmark
	./virtual_function | tee -a $(name).benchmark
	./pairlists | tee -a $(name).benchmark
