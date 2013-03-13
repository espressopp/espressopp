// ESPP_CLASS
#ifndef _VERSION_HPP
#define _VERSION_HPP
#include <string>
#include "boost/version.hpp"

#define MAJORVERSION 1
#define MINORVERSION 5
#define PATCHLEVEL   1
#include "hgversion.hpp"

namespace espresso {

  class Version {
  public:
    Version();
    std::string info();
    static void registerPython();

  private:
    int major;
    int minor;
    int patchlevel;
    std::string name;
    std::string hgrevision;
    std::string boostversion;
    std::string date;
    std::string time;
  };
}
#endif
