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
#include <stdlib.h>

// by the following define we make sure that logging with level << ERROR will
// be removed already at compile time; it has to be set before 
// we include the log4espp.hpp file 

#define LOG4ESPP_LEVEL_ERROR

// the next define makes sure that we overwrite a define that might have been
// set before or via -DLOG4ESPP_LEVEL_ERROR

#define LOG4ESPP_LEVEL_DEBUG

#include "logging/log4espp.hpp"

#define N 100

/**************************************************************************
*                                                                         *
*  MAIN program                                                           *
*                                                                         *
**************************************************************************/

LOG4ESPP_DEFINITION();  // in one unit put definitions of logging

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
   
        LOG4ESPP_DEBUG(rootLogger, i << " is not a prime number, can be divided by " << k);

     }

  }

  std::cout << "There are " << primeCounter << " primes up to " << N << "\n";

#ifdef LOG4ESPP_INFO_ENABLED
  if (LOG4ESPP_INFO_ON(rootLogger)) {

     int sum = 0;
     for (int k = 0; k < N; k++) sum += k;

     LOG4ESPP_INFO(rootLogger, "main program terminates with sum = " << sum);
  }
#endif


  if (LOG4ESPP_TRACE_ON(rootLogger)) {
    std::cout << "log4espp: trace enabled" << std::endl;
  }
  if (LOG4ESPP_DEBUG_ON(rootLogger)) {
    std::cout << "log4espp: debug enabled" << std::endl;
  }
  if (LOG4ESPP_INFO_ON(rootLogger)) {
    std::cout << "log4espp: info enabled" << std::endl;
  }
  if (LOG4ESPP_WARN_ON(rootLogger)) {
    std::cout << "log4espp: warn enabled" << std::endl;
  }
  if (LOG4ESPP_ERROR_ON(rootLogger)) {
    std::cout << "log4espp: error enabled" << std::endl;
  }
  if (LOG4ESPP_FATAL_ON(rootLogger)) {
    std::cout << "log4espp: fatal enabled" << std::endl;
  }

  LOG4ESPP_DEBUG(rootLogger, "log4espp: debug enabled");
  LOG4ESPP_TRACE(rootLogger, "log4espp: trace enabled");
  LOG4ESPP_INFO(rootLogger, "log4espp: info enabled");
  LOG4ESPP_WARN(rootLogger, "log4espp: warn enabled");
  LOG4ESPP_ERROR(rootLogger, "log4espp: error enabled");
  LOG4ESPP_FATAL(rootLogger, "log4espp: fatal enabled");
}
