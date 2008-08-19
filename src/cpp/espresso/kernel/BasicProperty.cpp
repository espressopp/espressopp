#include "BasicProperty.hpp"
#include "ParticleContainer.hpp"

#include <string>
#include <cassert>
#include <iostream>

using namespace std;

template<class X, class A> inline void Assert(A assertion) {
   if (!assertion) throw X();
}

BasicProperty::BasicProperty (int tsize, int nelems, const string& name) {

    typeSize = tsize;
    elemSize = nelems;
    propName = name;

    position  = -1;
    container = 0;
}

void BasicProperty::registerContainer(class ParticleContainer *pContainer, int pos) {

    Assert<std::bad_exception>(container == 0);

    // assert(container==0); 

    container = pContainer;
    position  = pos;

#ifdef DEBUG
    cout << "property " << propName << " registered for container " << pContainer << " at pos " << pos << "\n";
#endif

}

bool BasicProperty::isRegistered(class ParticleContainer *pContainer) {

#ifdef DEBUG
    cout << "property " << propName << " of container " << container << " check " << pContainer << "\n";
#endif

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

    cout << "BasicProperty (Name=" << propName << 
           ", typeSize=" << typeSize << ", size=" << (typeSize * elemSize) 
           << ", pos=" << position << ")\n"; 

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

