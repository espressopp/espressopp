#include "Simulation.hpp"

#include <string>
#include <cassert>
#include <iostream>

using namespace std;
using namespace boost;

Group::Group(string name) {

    this->name = name;
}

Group::~Group() {

    cout << "Group " << name << " will be deleted" << std::endl;
}

string Group::getName() {

    return name;
}

void Simulation::addGroup(GroupPtr g) {

    myGroups.push_back(g);
    cout << "Simulated: added group " << g->getName() << std::endl;
}

Simulation::Simulation() {
}

GroupPtr Simulation::getGroup(int n) {

   return myGroups[n];
}

Simulation::~Simulation() {

  cout << "Delete simulation" << std::endl;

    for (int i = 0; i < myGroups.size(); i++) {

        myGroups[i].reset();
    }

    cout << "Simulation has groups deleted, now free" << std::endl;
}

int Simulation::numberGroups() {

    return myGroups.size();

}

