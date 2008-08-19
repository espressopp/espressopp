#include "ParticleContainer.hpp"

#define DEFAULT_BLOCK_SIZE 256

#include <iostream>

/** \file ParticleContainer.cpp    i

Implementation of class ParticleContainer for storing particle data.

<b>Responsible:</b>
<a href="mailto:brandes@scai.fraunhofer.de">Thomas Brandes</a>

*/

/** This is the constructor of a ParticleContainer.
 *
 * \param numberBlocks specifies the initial number of blocks of particles
 * \param blockDefaultSize is the initial size used for a particle block
*/


/* dummy data to store just data for one particle / one property */

static int mockData [256];

ParticleContainer::ParticleContainer(int numberBlocks, int blockDefaultSize) {

}

void ParticleContainer::addProperty (BasicProperty *newProperty) {

}

void ParticleContainer::removeProperty (BasicProperty *p) {

}

ParticleRef ParticleContainer::addParticle (int blockIndex, int key) {

  return ParticleRef(0,this);
}

void* ParticleContainer::getAllParticleData(BasicProperty *prop) {

  return NULL;
}

void* ParticleContainer::getParticleDataPtr(ParticleRef pRef, BasicProperty *prop) {

  return mockData;
}

void ParticleContainer::setParticleData (ParticleRef pRef, BasicProperty *prop, void *data) {

  int size = prop->getSize();
  memcpy (mockData, data, size);
  
}

void ParticleContainer::getParticleData(ParticleRef pRef, BasicProperty *prop, void *data) {

  int size = prop->getSize();
  memcpy (data, mockData, size);
  
}

void ParticleContainer:: moveParticle(ParticleRef p, int blockIndex) {

}

int ParticleContainer::getParticleKey(ParticleRef pRef) {

}

int* ParticleContainer::getAllParticleKey() {

}

ParticleRef ParticleContainer::findParticle (int key) {

  return ParticleRef(0,this);

}

void ParticleContainer::printInfo () {

   cout << "--begin-- ParticleContainer::printInfo\n";
   cout << "---end--- ParticleContainer::printInfo\n";
}
