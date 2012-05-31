// ESPP_CLASS
#ifndef _INTEGRATOR_TDFORCE_HPP
#define _INTEGRATOR_TDFORCE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Real3D.hpp"
#include "SystemAccess.hpp"
#include "interaction/Interpolation.hpp"
#include <map>


#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"


namespace espresso {
  namespace integrator {

    class TDforce : public Extension {

      public:
        TDforce(shared_ptr<System> system);
        ~TDforce();

        /** Setter for the filename, will read in the table. */
        void addForce(int itype, const char* _filename, int type);
        const char* getFilename() const { return filename.c_str(); }

        void applyForce();

        void setCenter(real x, real y, real z);

        static void registerPython();

      private:

        boost::signals2::connection _applyForce;

        void connect();
        void disconnect();

        Real3D center;
        std::string filename;
        typedef shared_ptr <interaction::Interpolation> Table;
        std::map<int, Table> forces; // map type to force

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
