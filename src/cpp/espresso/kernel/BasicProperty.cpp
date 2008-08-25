#include "BasicProperty.hpp"
#include "ParticleContainer.hpp"

#include <string>
#include <cassert>
#include <iostream>

using namespace std;

// Define a static logger for this class

LOG4ESPP_LOGGER(BasicProperty::logger, "kernel.BasicProperty");

template<class X, class A> inline void Assert(A assertion) {
   if (!assertion) throw X();
}

BasicProperty::BasicProperty (int _typeSize, int _elemSize, const string& _propName) {

    typeSize = _typeSize;
    elemSize = _elemSize;
    propName = _propName;

    position  = -1;
    container = 0;

    LOG4ESPP_INFO(logger,"property " << propName << " defined, size = " << typeSize * elemSize);
}

void BasicProperty::registerContainer(class ParticleContainer *pContainer, int pos) {

    Assert<std::bad_exception>(container == 0);

    container = pContainer;
    position  = pos;

    LOG4ESPP_INFO(logger, "property " << propName << " registered for container " 
                           << pContainer << " at pos " << pos);

}

bool BasicProperty::isRegistered(class ParticleContainer *pContainer) {

    LOG4ESPP_DEBUG(logger,"property " << propName << " of container " << container << " check " << pContainer);

    return (pContainer == container);
}

void BasicProperty::updatePosition(int pos) {

    position = pos;
}

int BasicProperty::getPosition() {

    return position;
}

int BasicProperty::getSize() {

    return typeSize*elemSize;
}

string BasicProperty::getName() {

    return propName;
}

void BasicProperty::deregister() {

    assert(container!=0); 
    container = 0;
}

void BasicProperty::printInfo() {

    LOG4ESPP_INFO(logger, "BasicProperty (Name=" << propName << 
                          ", typeSize=" << typeSize << ", size=" << (typeSize * elemSize) 
                          << ", pos=" << position); 

}

void BasicProperty::printData(void *data) {

    char *ptr = (char *) data;

    cout << "BasicProperty (Name=" << propName << 
           ", typeSize=" << typeSize << ", size=" << (typeSize * elemSize) 
           << ", pos=" << position << ")\n"; 

    cout << "data: ";

    for (int i = 0; i < typeSize * elemSize; i++) {
       cout << (int) *ptr++  << " ";
    }
 
    cout << "\n";
}

BasicProperty::~BasicProperty() {
	
}

void BasicProperty::storeData(ParticleRef pRef, void *data) {

    assert(container!=0); 
    
    container->setParticleData (pRef, this, data);
}

void BasicProperty::fetchData(ParticleRef pRef, void *data) {

    assert(container!=0); 
    
    container->getParticleData (pRef, this, data);
}

void * BasicProperty::getPtrData(ParticleRef pRef) {

    assert (container!=0);
 
    container->getParticleDataPtr (pRef, this);
}

