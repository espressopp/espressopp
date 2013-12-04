// ESPP_CLASS
#ifndef _ADRESS_HPP
#define _ADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espresso {

  namespace integrator {

      class Adress : public Extension {

      public:
        shared_ptr<VerletListAdress> verletList;
        shared_ptr<FixedTupleListAdress> fixedtupleList;
        bool KTI;
        
        real dhy;
        real pidhy2;
        real dex;
        real dex2;
        real dexdhy;
        real dexdhy2;
        
        Adress(shared_ptr<System> _system, shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, bool _KTI = false);

        ~Adress();
        
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        boost::signals2::connection _SetPosVel, _initForces, _integrate1, _integrate2, _aftCalcF;
        
        void integrate1(real&);
        void initForces();
        void SetPosVel();
        void integrate2();
        void aftCalcF();

        void connect();
        void disconnect();
        
        real weight(real);
        real weightderivative(real);
      };

  }

}

#endif
