/**************************************************************************
*                                                                         *
*  Author      : Dr. Thomas Brandes, SCAI, FhG                            *
*  Copyright   : SCAI, FhG St. Augustin, Germany                          *
*  Date        : Jul 08                                                   *
*  Last Update : Jul 08                                                   *
*                                                                         *
*  Sample main program of using Espresso++                                *
*                                                                         *
*  Module      : main                                                     *
*                                                                         *
*  Function    : int main (int argc, char **argv)                         *
*                                                                         *
**************************************************************************/

#include <string>
#include <iostream>
#include <stdlib.h>

#include "BasicProperty.hpp"
#include "ParticleContainer.hpp"
#include "log4espp.hpp"

int main (int argc, char **argv) {

  LOG4ESPP_CONFIGURE();       // read runtime configuration for logging

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  int bigBlocks = 1000;
  int sizeBigBlocks = 10000;

  ParticleContainer *container = new ParticleContainer(bigBlocks, sizeBigBlocks);

  ArrayProperty<double,3> acceleration ("acceleration");
  ArrayProperty<double,3> position     ("position");

    // for registration of a property we pass it by reference

  container->addProperty (&acceleration);
  container->addProperty (&position);

  int key = 1131;

  LOG4ESPP_INFO(rootLogger, "define " << bigBlocks << " blocks with " 
                                      << sizeBigBlocks << " particles each");

  for (int block = 0; block < bigBlocks; block++) {

    if ((block & 15) == 0) {
      LOG4ESPP_DEBUG(rootLogger, "define block " <<  block);
    }

    for (int j = 0; j < sizeBigBlocks; j++) {

      double val [3];
      val [0] = 1.0;
      val [1] = -1.0;
      val [2] = 0.0;
      ParticleRef refParticle = container->addParticle (block, key++);
      acceleration.setData (refParticle, val);
      
      double (&pos) [3] = position [refParticle];

      pos [0] = drand48();
      pos [1] = drand48();
      pos [2] = drand48();
   }
  }

  LOG4ESPP_INFO(rootLogger, "All particles inserted");

   // loop over all particles sum up accelerations

   double sumacc [3];
   double sumpos [3];

   sumacc[0] = 0.0;
   sumacc[1] = 0.0;
   sumacc[2] = 0.0;
   sumpos[0] = 0.0;
   sumpos[1] = 0.0;
   sumpos[2] = 0.0;

   int count = 0;

   LOG4ESPP_INFO(rootLogger, "now loop over all particles");

   for (ParticleIterator it = container->begin(); it !=container->end(); ++it) {

      ParticleRef refParticle = *it;

      sumacc[0] += acceleration[refParticle][0];
      sumacc[1] += acceleration[refParticle][1];
      sumacc[2] += acceleration[refParticle][2];
      sumpos[0] += position[refParticle][0];
      sumpos[1] += position[refParticle][1];
      sumpos[2] += position[refParticle][2];

      count++;
   }

   LOG4ESPP_INFO(rootLogger,"looped over " << count << " particles");
   LOG4ESPP_INFO(rootLogger,"sum acceleration = " << sumacc[0] << ", " << sumacc[1] << ", " << sumacc[2]);
   LOG4ESPP_INFO(rootLogger,"sum positions    = " << sumpos[0] << ", " << sumpos[1] << ", " << sumpos[2]);

}
