#ifndef BasicProperty_H
#define BasicProperty_H

/** \file BasicProperty.hpp    Class for basic properties of particles.

<b>Responsible:</b>
<a href="mailto:brandes@scai.fraunhofer.de">Thomas Brandes</a>

A property is defined by its name, by its type and the number of elements of this type.

A property can belong to exactly one property list. In this case the property gets 
a position that corresponds to the position in the list.

*/

#include <string>
#include <iostream>

class ParticleRef {

public:

  int   index;                  // index position of particle
  class ParticleContainer *pc;  // particle container itself

  ParticleRef (int index, ParticleContainer *pc) {

     this->index = index;
     this->pc    = pc;
  }
};

class BasicProperty {

protected:

    int    elemSize;  /* number of elements       */

    std::string propName;  /* name of the property     */

private:

    int    typeSize;  /* type size                */

    class  ParticleContainer *container;  /* reference to the particle container to which property belongs */

    int    position;  /* row position in ParticleContainer */

public:

/**
   Constructor of a basic property.

   @param tsize specifies the number of bytes for the used type, e.g. 4 for int.
   @param nelem specifies the number of elements used in the property.
   @param propName is a string that specifies the name of the property.
   @return a new basic property.
*/
 
    BasicProperty(int tsize, int nelem, const std::string& propName);

/**
   Ask the property how many bytes are needed for a value of this property.

   @return number of bytes as int value.
 */

    int getSize(); 

/**
   Ask for the name of the property.

   @return a string with the name of this property.
 */

    std::string getName();

/**
   Register this property for a particle containter. By this way the
   property can be used to store and fetch data.

   @param pContainer is a pointer to the particle container.
   @param pos specifies the position of the attribute in the container.
 */

    void registerContainer(class ParticleContainer *pContainer, int pos);

/**
   Predicate to prove whether this property is really registered for
   the given container. 

   @param pContainer is a pointer to the particle container.
   @return true if this property is registered for pContainer
 */

    bool isRegistered(class ParticleContainer *pContainer);

    void updatePosition(int pos);

    int  getPosition();

    void deregister();

    void printInfo();

    void storeData(ParticleRef pRef, void *data);

    void fetchData(ParticleRef pRef, void *data);

    void* getPtrData (ParticleRef pRef);

    virtual void printData(void *data);

    virtual ~BasicProperty();
    
};

/**************************************************************************
*                                                                         *
*  class ScalarProperty<T>                                                *
*                                                                         *
**************************************************************************/

template <class T> 
class ScalarProperty : public BasicProperty {

public :

   T myDefaultValue;

/**
   Constructor for a scalar property.

   @param propName is a string to identify the property.
   @param inputVal will be the default value for this property.
 */

   ScalarProperty (std::string propName, T inputVal);

/**
   Constructor for a scalar property.

   @param propName is a string to identify the property.
 */

   ScalarProperty (std::string propName);

/**
   This routine sets for a particle the default value.

   @param ref specifies the particle for which we set the property.
 */

   void setDefault(ParticleRef ref);

/**
   Setter routine for a scalar property. The property must
   have been registered for a particle container.

   @param ref specifies the particle that will be referenced.
   @param val is the value of the property to be set.
 */

   void setData(ParticleRef ref, T val);

/**
   Getter routine for a scalar property. The property must
   have been registered for a particle container.

   @param ref specifies the particle that will be referenced.
   @return scalar value of the property for the specified particle.
 */

   T getData(ParticleRef ref);

/**
   Index operator to get for a particle reference the
   reference to the corresponding type.

   @param ref specifies the particle that will be referenced.
   @return Reference to the property data of the particle.
 */

   T&  operator[] (ParticleRef ref);

};

template<class T>
ScalarProperty<T>::ScalarProperty (std::string propName)

   : BasicProperty(sizeof(T), 1, propName)

{  myDefaultValue = 0;
}

template<class T>
ScalarProperty<T>::ScalarProperty (std::string propName, T inputVal)

   : BasicProperty(sizeof(T), 1, propName)

{  myDefaultValue = inputVal;
}

template<class T>
void ScalarProperty<T>::setData (ParticleRef ref, T val) {

    storeData(ref, (void *) &val);
}

template<class T>
void ScalarProperty<T>::setDefault (ParticleRef ref) {

    storeData(ref, (void *) &myDefaultValue);
}

template<class T>
T ScalarProperty<T>::getData (ParticleRef ref) {

    T result;
    fetchData(ref, (void *) &result);
    return result;
}

template<class T>
T& ScalarProperty<T>::operator[] (ParticleRef ref) {

    T* data = (T*) getPtrData(ref);
    return  data[0];
}

/**************************************************************************
*                                                                         *
*  class ArrayProperty<T>                                                 *
*                                                                         *
**************************************************************************/

template <class T, int n> 
class ArrayProperty : public BasicProperty {

public :

    ArrayProperty (std::string propName);

    void setDefault (ParticleRef ref);

    void setData(ParticleRef ref, T inputVal[n]);

    void getData(ParticleRef ref, T outputVal[n]);

    T (& operator[] (ParticleRef ref))[n];
};

template<class T, int n>
ArrayProperty<T,n>::ArrayProperty (std::string propName)

   : BasicProperty(sizeof(T), n, propName)

{  
}

template<class T, int n>
void ArrayProperty<T,n>::setData (ParticleRef ref, T inputVal[n]) {

   storeData(ref, inputVal);
}

template<class T, int n>
void ArrayProperty<T,n>::setDefault (ParticleRef ref) {

   T val[n];

   for (int i = 0; i < n; i++) {
       val [i] = 0;
   }

   storeData(ref, val);
}

template<class T, int n>
void ArrayProperty<T,n>::getData (ParticleRef ref, T outputVal[n]) {

   fetchData(ref, outputVal);
}

template<class T, int n>
T (& ArrayProperty<T,n>::operator[] (ParticleRef ref))[n] {

   // illegal cast: T (& data) [n] = (T (&) [n]) getPtrData(ref);

   T (* data) [n] = (T (*) [n]) getPtrData(ref);

   return data[0];
}

#endif
