#include <boost/python.hpp>
#include <boost/shared_ptr.hpp> 

#include <iostream>
#include "Simulation.hpp"

using namespace boost::python;
using namespace boost;

BOOST_PYTHON_MODULE(_Simulation)

{    class_<Simulation> ("Simulation", init<>())
    .def("addGroup", &Simulation::addGroup)
    .def("getGroup", &Simulation::getGroup)
    .def("numberGroups", &Simulation::numberGroups)
    ;

    class_<Group, GroupPtr >("Group", init<const char *>())
    .def("getName", &Group::getName)
    ;

}

