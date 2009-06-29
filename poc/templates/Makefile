CXXFLAGS=-O3 -g

all: multi_parameter_template

.PHONY: clean all FORCE

FORCE:

multi_parameter_template: multi_parameter_template.C classes.C instantiations.C FORCE
	time -p $(CXX) $(CXXFLAGS) $< 

classes.C: generate_classes.pl
	perl generate_classes.pl > classes.C

instantiations.C: generate_instantiations.pl
	perl generate_instantiations.pl > instantiations.C

clean:
	rm -f classes.C instantations.C 
	rm -f multi_parameter_template
