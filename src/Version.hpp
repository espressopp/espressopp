// ESPP_CLASS
#ifndef _VERSION_HPP
#define _VERSION_HPP
#include <string>

#define MAJORVERSION 0
#define MINORVERSION 9
#define PATCHLEVEL   0
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
    std::string date;
    std::string time;
  };
}
#endif
