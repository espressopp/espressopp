#ifndef Simulation_H
#define Simulation_H

#include <string>
#include <vector>

#include "boost/smart_ptr.hpp" 

class Group;

typedef boost::shared_ptr<Group> GroupPtr;

class Group {

private:

    std::string name;
 
public:

    Group(std::string name);
    ~Group();
    std::string getName();
};

class Simulation {

private:
 
    std::vector<GroupPtr> myGroups;

public:
 
    Simulation();
    ~Simulation();

    void addGroup(GroupPtr g);
    GroupPtr getGroup(int);
    int numberGroups();

};

#endif
