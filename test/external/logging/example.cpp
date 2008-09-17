/**************************************************************************
*                                                                         *
*  Author      : Dr. Thomas Brandes, SCAI, FhG                            *
*  Copyright   : SCAI, FhG St. Augustin, Germany                          *
*  Date        : Sep 08                                                   *
*  Last Update : Sep 08                                                   *
*                                                                         *
*  Sample program to demonstrate LOGGING facilities                       *
*                                                                         *
*  Module      : example                                                  *
*                                                                         *
**************************************************************************/

#include <iostream>

#include "log4espp.hpp"

#define N 100

/**************************************************************************
*                                                                         *
*  MAIN program                                                           *
*                                                                         *
**************************************************************************/

int main (int argc, char **argv) {

  int primeCounter = 0;

  int i;

  LOG4ESPP_CONFIGURE();       // read runtime configuration for logging

  LOG4ESPP_ROOTLOGGER(rootLogger);  // get the rootLogger

  LOG4ESPP_INFO(rootLogger, "main program starts");

  for (i = 2; i <= N; i++) {

      int is_prime = 1;
      int k = 2; 

      while (k*k <= i) {

        if (i % k == 0) {

          is_prime = 0;
          break;
        }

        k++;
 
      }

     if (is_prime) {

        primeCounter ++;

        LOG4ESPP_DEBUG(rootLogger, i << " is a prime number");

     } else {
   
        LOG4ESPP_DEBUG(rootLogger, i << " is not a prime number");

     }

  }

  std::cout << "There are " << primeCounter << " primes up to " << N << "\n";

  LOG4ESPP_INFO(rootLogger, "main program terminates");

}
