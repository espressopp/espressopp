/**************************************************************************
*                                                                         *
*  Author      : Dr. Thomas Brandes, SCAI, FhG                            *
*  Copyright   : SCAI, FhG St. Augustin, Germany                          *
*  Date        : Jul 08                                                   *
*  Last Update : Jul 08                                                   *
*                                                                         *
*  Sample program to calculate the circle of Halley's comet               *
*                                                                         *
*  Module      : halley                                                   *
*                                                                         *
*  Function    : int halley (int argc, char **argv)                       *
*                                                                         *
**************************************************************************/

#include <string>
#include <iostream>
#include <cmath>

#include "BasicProperty.hpp"
#include "ParticleContainer.hpp"
#include "log4espp.hpp"

/**************************************************************************
*                                                                         *
*  Global definitions                                                     *
*                                                                         *
**************************************************************************/

#define PLANET_SUN     0
#define PLANET_EARTH   1
#define PLANET_JUPITER 2
#define PLANET_HALLEY  3

#define sqr(x) ((x) *(x))

#define REAL double
#define DIM  2

/**************************************************************************
*                                                                         *
*  compForce                                                              *
*                                                                         *
**************************************************************************/

void compForce(ParticleContainer *p,
               ArrayProperty<REAL,DIM> *F,
               ScalarProperty<REAL> *m,
               ArrayProperty<REAL,DIM> *x)   {

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  for (ParticleIterator it1 = p->begin(); it1 !=p->end(); ++it1) {
   
    ParticleRef p1 = *it1;

    REAL &m1     = (*m)[p1];
    REAL(&F1)[DIM] = (*F)[p1];
    REAL(&x1)[DIM] = (*x)[p1];

    for (int i = 0; i < DIM; i++) {
      F1[i] = 0.0;
    }

    for (ParticleIterator it2 = p->begin(); it2 !=p->end(); ++it2) {

       ParticleRef p2 = *it2;

       if (it1 != it2) {

         REAL  &m2       = (*m)[p2];
         REAL (&F2)[DIM] = (*F)[p2];
         REAL (&x2)[DIM] = (*x)[p2];

         // r for the square root distance

         REAL r = 0.0;
         for (int i = 0; i < DIM; i++) {
            r += sqr(x1[i]- x2[i]);
         }
 
         REAL f = m1 * m2 / (sqrt(r) * r);

         /* 
         LOG4ESPP_DEBUG(rootLogger, "Force f = " << f << 
                        " r = " << r << ", m1 = " << m1 << ", m2 = " << m2);
         */

         for (int i = 0; i < DIM; i++) {
            F1[i] += f * (x2[i] - x1[i]);
         }
       }
    }
  }

  for (ParticleIterator it = p->begin(); it !=p->end(); ++it) {

    ParticleRef pp = *it;

    REAL(&pF)[DIM]    = (*F)[pp];

    LOG4ESPP_DEBUG(rootLogger, "Force " << p->getParticleKey(pp) 
                   << " = " << pF[0] << ", " << pF[1]);

  }
}

void compPosition(ParticleContainer *pC, REAL deltaTime,
                  ArrayProperty<REAL,DIM> *F,
                  ArrayProperty<REAL,DIM> *Fold,
                  ScalarProperty<REAL> *m,
                  ArrayProperty<REAL,DIM> *x,
                  ArrayProperty<REAL,DIM> *v) {

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  for (ParticleIterator it = pC->begin(); it !=pC->end(); ++it) {

    ParticleRef refParticle = *it;

    REAL(&px)[DIM]    = (*x)[refParticle];
    REAL(&pv)[DIM]    = (*v)[refParticle];
    REAL&pm           = (*m)[refParticle];
    REAL(&pF)[DIM]    = (*F)[refParticle];
    REAL(&pFold)[DIM] = (*Fold)[refParticle];

    REAL a = deltaTime * 0.5 / pm;
    
    for (int i = 0; i < DIM; i++) {
      px[i] += deltaTime * (pv[i] + a * pF[i]);
      pFold[i] = pF[i];
    }
   
    LOG4ESPP_DEBUG(rootLogger, "pos " << pC->getParticleKey(refParticle) 
                          << " = " << px[0] << ", " << px[1]);
  }
}

void compVelocity(ParticleContainer *p, 
                  REAL deltaTime,
                  ArrayProperty<REAL,DIM> *F,
                  ArrayProperty<REAL,DIM> *Fold,
                  ScalarProperty<REAL>    *m,
                  ArrayProperty<REAL,DIM> *v) {

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  for (ParticleIterator it = p->begin(); it !=p->end(); ++it) {

    ParticleRef refParticle = *it;

    REAL(&pv)[DIM]    = (*v)[refParticle];
    REAL&pm           = (*m)[refParticle];
    REAL(&pF)[DIM]    = (*F)[refParticle];
    REAL(&pFold)[DIM] = (*Fold)[refParticle];

    REAL a = deltaTime * 0.5 / pm;

    for (int i = 0; i < DIM; i++) {
       pv[i] += a * (pF[i] + pFold[i]);
    }

    LOG4ESPP_DEBUG(rootLogger, "velocity " << p->getParticleKey(refParticle) 
                   << " = " << pv[0] << ", " << pv[1]);
  }

}

void compStatistic(ParticleContainer *pC,
                   REAL time,
                   ScalarProperty<REAL> *m,
                   ArrayProperty<REAL,DIM> *v) {

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  REAL e = 0.0;

  for (ParticleIterator it = pC->begin(); it !=pC->end(); ++it) {

      ParticleRef refParticle = *it;
      REAL(&pv)[DIM] = (*v)[refParticle];
      REAL v = 0.0;
      for (int i = 0; i < DIM; i++) {
         v += sqr(pv[i]);
      }
      e += 0.5 * (*m)[refParticle] * v;
  }
   
  LOG4ESPP_INFO(rootLogger, "kinetic energy at time " << time << " is " << e);
}

void outputResults(ParticleContainer *pC,
                   ArrayProperty<REAL,DIM> *x) {

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  for (ParticleIterator it = pC->begin(); it !=pC->end(); ++it) {
      ParticleRef refParticle = *it;
      REAL(&px)[DIM] = (*x)[refParticle];
      LOG4ESPP_INFO(rootLogger, "pos of " << pC->getParticleKey(refParticle) 
                    << " = " << px[0] << ", " << px[1]);
  }
}

/**************************************************************************
*                                                                         *
*  Time integration                                                       *
*                                                                         *
**************************************************************************/

void timeIntegration (REAL startTime, REAL deltaTime, 
                      REAL endTime, ParticleContainer *p,
                      ArrayProperty<REAL,DIM> *F,
                      ArrayProperty<REAL,DIM> *Fold,
                      ScalarProperty<REAL> *m,
                      ArrayProperty<REAL,DIM> *x,
                      ArrayProperty<REAL,DIM> *v)

{ REAL time = startTime;

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  compForce (p, F, m, x);

  int itCounter = 0;

  REAL printDeltaTime = 10.0;

  REAL nextPrintTime = startTime + printDeltaTime;

  while (time < endTime) {
   
      time += deltaTime;

      compPosition(p, deltaTime, F, Fold, m, x, v);
      compForce(p, F, m, x);
      compVelocity(p, deltaTime, F, Fold, m, v);
      
      if (time >= nextPrintTime) {

         nextPrintTime += printDeltaTime;
         compStatistic(p, time, m, v);
         outputResults(p, x);
      }

      itCounter ++;
  }
}
      
int main (int argc, char **argv) {

  LOG4ESPP_CONFIGURE();       // read runtime configuration for logging

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  // one block with, 4 planets 

  ParticleContainer *container = new ParticleContainer(1, 4);

  ArrayProperty<REAL,DIM> F ("force");
  ArrayProperty<REAL,DIM> Fold ("oldforce");
  ScalarProperty<REAL> m ("mass");
  ArrayProperty<REAL,DIM> x ("position");
  ArrayProperty<REAL,DIM> v ("velocity");

  // for registration of a property we pass it by reference

  container->addProperty (&F);
  container->addProperty (&Fold);
  container->addProperty (&m);
  container->addProperty (&x);
  container->addProperty (&v);

  int block = 0;  // all planets belong to the first block

  ParticleRef sun = container->addParticle (block, PLANET_SUN);
  ParticleRef earth = container->addParticle (block, PLANET_EARTH);
  ParticleRef jupiter = container->addParticle (block, PLANET_JUPITER);
  ParticleRef halley  = container->addParticle (block, PLANET_HALLEY);

  // set masses of the planets

  m[sun]     = 1.0;
  m[earth]   = 3.0e-6;
  m[jupiter] = 9.55e-4;
  m[halley]  = 1e-14;

  // set positions of the planets

  x[sun][0] = 0.0;
  x[sun][1] = 0.0;
  x[earth][0] = 0.0;
  x[earth][1] = 1.0;
  x[jupiter][0] = 0.0;
  x[jupiter][1] = 5.36;
  x[halley][0] = 34.75;
  x[halley][1] = 0.0;

  // set velocity of the planets

  v[sun][0] = 0.0;
  v[sun][1] = 0.0;
  v[earth][0] = -1.0;
  v[earth][1] = 0.0;
  v[jupiter][0] = -0.425;
  v[jupiter][1] = 0.0;
  v[halley][0] = 0.0;
  v[halley][1] = 0.0296;

  LOG4ESPP_INFO(rootLogger, "All planets inserted");

  REAL startTime = 0.0;
  REAL deltaTime = 0.015; 
  REAL endTime   = 468.5; 

  LOG4ESPP_INFO(rootLogger, "start time integration");

  timeIntegration (startTime, deltaTime, endTime, container, &F, &Fold, &m, &x, &v);

  LOG4ESPP_INFO(rootLogger, "end time integration");

  LOG4ESPP_INFO(rootLogger, "Finished");
}
