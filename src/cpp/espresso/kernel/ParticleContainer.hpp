#ifndef ParticleContainer_H
#define ParticleContainer_H

#include <string>
#include <vector>

#include "BasicProperty.hpp"
#include "logging.hpp"

/** \file ParticleContainer.hpp    Class for storing particle data.

<b>Responsible:</b>
<a href="mailto:brandes@scai.fraunhofer.de">Thomas Brandes</a>

A particle container stores the data of its particles for the registered properties
(e.g. position, force, charge). 

The particles are organized in blocks where a block can be considered
as a dynamically growing array of particles. There is no order of particles within a block.
A block of particles might be used to represent the particles on one cell 
(domain decomposition) or to hold the data of particles from 
another processor (atomic decomposition).

From the users point of view, the user can be sure that for each property the data 
of all particles within one block is stored contiguously in the memory.
Therefore, moving a particle from one block into another block implies moving of all
its property data in the memory.

A particle container itself has no knowledge about any distribution of particles among the
worker. Each worker has its own particle container. The information about the 
organization of the blocks is left to the use of this class.

*/

using namespace std;

#define BLOCK_DEFAULT_SIZE 256

/** \class ParticleIterator  

    Iterator class for traversing all particles in a ParticleContainer.

*/

class ParticleIterator {

private:

   ParticleContainer *myContainer;  /** reference to the corresponding container */
   int  block;                      /** actual block index                       */
   int  blockpos;                   /** actual position in the block             */
   int  pref;                       /** actual particle index                    */

public:

   friend class ParticleContainer;

/** This routine sets the iterator to the next valid position of a particle.
*/
   void nextPos();

   bool operator== (const ParticleIterator& i) const 

      { return (block == i.block) && (blockpos == i.blockpos); }

   bool operator!= (const ParticleIterator& i) const

      { return (block != i.block) || (blockpos != i.blockpos); }

   ParticleRef operator* ();           // Gives a ref to an actual element

   ParticleIterator& operator++ ();    // Moves forward

   // ParticleIterator  operator++ (int);     // Moves forward
   // ParticleIterator& operator-- ();        // Moves backwards
   // ParticleIterator  operator-- (int);     // Moves backwards

};

/** \class ParticleContainer

    An object of this class can be used to store particle data for a
    certain number of properties. The particles are organized in blocks.

*/

class ParticleContainer {

private:

    int noBlocks;              // number of blocks in the container

    vector<int> blockOffset;   // offset for each block in the container data

    vector<int> blockSize;     // current number of particles in the block

    vector <BasicProperty*> propertyList;  // list of all basic properties 

    // Hash table to map the keys to indexes

    vector <void*> propertyData;

    int *keyData;              // pointer to the memory used for the key values

    int allocatedSize; 

    static LOG4ESPP_DECL_LOGGER(myLogger);

public:
 
/** This constructor creates a particle container for a number of blocks with a default size.
    The data will be allocated as soon as properties are added.

 * \param numberBlocks specifies the initial number of blocks of particles
 * \param blockDefaultSize is the initial size used for a particle block
*/


    ParticleContainer(int numberBlocks, int blockDefaultSize = BLOCK_DEFAULT_SIZE);

    // add a property to particles, particle container will allocate
    // data for this property

    void addProperty (BasicProperty *newProperty);

    // remove a property to particles, particle container will free
    // data for this property

    void removeProperty (BasicProperty *p);

    // add a particle for a block, 0 <= blockIndex < numberBlocks 
    // if one block does not contain enough data 
    // IMPORTANT: there is never any assumption about the order in the  block

/** This routine is used to add a completely new particle in the particle container.

 * \param blockIndex specifies the block in which the particle is created
 * \param key must be a unique integer that is used to identify the particle.

*/

    ParticleRef addParticle (int blockIndex, int key);
    
    // update particle data in the particle container 
    // where the value is taken from p->myVal
    
    void setParticleData (ParticleRef pRef, BasicProperty *p, void *data);

    // ask for a special property of a particle and return it as prop->myVal 

    void getParticleData (ParticleRef pRef, BasicProperty *p, void *data);

    // ask for the system key of the particle

    int getParticleKey (ParticleRef pRef);

    // move one particle in another block (key remains)

    void moveParticle (ParticleRef p, int blockIndex);

    // Pick up an array containg the property data of all data; 
    // this array can be addressed with the particle Reference

    void* getAllParticleData (BasicProperty *p);

    // get pointer to the data of one particle

    void* getParticleDataPtr (ParticleRef pRef, BasicProperty *p);

    // get an integer array containing all keys 

    int* getAllParticleKey ();

    // Pick up an integer array containing all the keys of the particles
    // this array can be addressed with the particle Reference

    // find a particle by a key; this routine is assumed to be fast by hashing.

    ParticleRef findParticle (int key);

    // print info about this container in stdout

    void printInfo();

    ParticleIterator begin();   // Pointer to the first particle
    ParticleIterator end();     // Pointer to the last particle

    friend class ParticleIterator;

};

#endif
