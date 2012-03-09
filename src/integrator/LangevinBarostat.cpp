#include "python.hpp"
#include "LangevinBarostat.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

#include "esutil/RNG.hpp"

#include "analysis/Pressure.hpp"
#include "bc/BC.hpp"

#include "mpi.hpp"

namespace espresso {

  using namespace std;
  using namespace analysis;
  
  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(LangevinBarostat::theLogger, "LangevinBarostat");

    LangevinBarostat::LangevinBarostat(shared_ptr<System> _system, shared_ptr< esutil::RNG > _rng) : SystemAccess(_system), rng(_rng){
      // external parameters
      gammaP = 0.0;
      mass = 0.0;
      externalPressure = 0.0;
      
      // internal parameter
      momentum = 0.0;

      //rng = _rng;
      
      LOG4ESPP_INFO(theLogger, "LangevinBarostat constructed");
    }

    void LangevinBarostat::setGammaP(real _gammaP){
      gammaP = _gammaP;
    }
    real LangevinBarostat::getGammaP(){
      return gammaP;
    }

    void LangevinBarostat::setPressure(real _press){
      externalPressure = _press;
    }
    real LangevinBarostat::getPressure(){
      return externalPressure;
    }
    
    void LangevinBarostat::setMass(real _mass){
      mass = _mass;
    }
    real LangevinBarostat::getMass(){
      return mass;
    }

    LangevinBarostat::~LangevinBarostat(){}

    void LangevinBarostat::updVolume(real dt_2){
      System& system = getSystemRef();
      
      // V(t+1/2*dt) = V(t) + 1/2*dt*V';    V' = d*V*pe/W
      real scale_factor = pow( 1 + dt_2 * 3.0 * momentum, 1./3.);  // calculating the current scaling parameter
      
      //mpi::communicator communic = *system.comm;
      //cout << "VOL:    " << 1 + dt_2 * 3.0 * momentum << "     cpu:"<< communic.rank() <<  "         !!!!!!!!!!!!!" << endl ;

      system.scaleVolume( scale_factor, false);
    }
    
    // @TODO should be optimized 
    void LangevinBarostat::updVolumeMomentum(real dt_2){
      // dt_2 is timestep/2. 
      System& system = getSystemRef();
      Real3D Li = system.bc -> getBoxL(); // getting the system size
      real V = Li[0] * Li[1] * Li[2];     // system volume
      
      // get a random value for each vector component
      mpi::communicator communic = *system.comm;
      real rannum; 
      if (communic.rank() == 0) rannum = rng->normal();
      boost::mpi::broadcast(communic, rannum, 0);
      
//*************************************************************************************************
      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      real v2sum;
      real v2 = 0.0;
      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        const Particle& p = *cit;
        v2 = v2 + p.mass() * (p.velocity() * p.velocity());
      }
      boost::mpi::all_reduce( communic, v2, v2sum, std::plus<real>());
//*************************************************************************************************
      
      static Pressure Pcurrent( getSystem() );
      real P = Pcurrent.compute();  // calculating the current pressure in system
      
      // @TODO one should check that X is not the instantaneous pressure, 
      // while it does not include the white noise from thermostat
      real X = P; 
      // local momentum derivative;   3 - dimensions
      real pe_deriv = 3 * V * (X - externalPressure)  + pref6 * v2sum +
              pref4 * momentum + pref5 * rannum;
      
      pe_deriv /= mass; // momentum is already normalized by mass
      
      momentum += dt_2 * pe_deriv;
    }
    
    void LangevinBarostat::updForces(){
      LOG4ESPP_DEBUG(theLogger, "barostating");

      System& system = getSystemRef();
      CellList cells = system.storage->getRealCells();

      for(CellListIterator cit(cells); !cit.isDone(); ++cit){
        frictionBarostat(*cit);
      }
    }
    
    real LangevinBarostat::updDisplacement(){
      // should return pe/W
      return (momentum); // momentum is normalized by fictitious mass already
    }

    void LangevinBarostat::frictionBarostat(Particle& p){
      
      p.force() += pref3  * momentum * p.velocity() * p.mass();
      
      LOG4ESPP_TRACE(theLogger, "new particle force = " << p.force());
    }

    void LangevinBarostat::initialize(real timestep, real desiredTemperature){
      // calculate the prefactors
      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
              ", gammaP = " << gammaP << 
              ", external pressure = " << externalPressure <<
              ", fictitious mass = " << mass );

      // barostat
      /* Nf - degrees of freedom. For N particles in d-dimentional system without constrains
       * Nf = d * N
       * !!!!!!!!!!!!!!!!!!!!! @TODO Should be changed, now just test case !!!!!!!!!!!
       */
      // determine number of local particles and total particles
      System& system = getSystemRef();
      int N, Nsum;
      N = system.storage->getNRealParticles();
      boost::mpi::all_reduce(*mpiWorld, N, Nsum, std::plus<int>());
      real Ntotal = Nsum;

      pref6 = 1./Ntotal;
      
      pref3 = - (1+pref6);
      
      pref4 = -gammaP;
      // normal distribution prefactor
      pref5 = sqrt(24.0 * desiredTemperature * gammaP * mass / timestep);

//      mpi::communicator communic = *system.comm;
//      cout << "g:  " << gammaP<< "   P: "<< externalPressure << "    M: " << mass<< "   T:  "<< desiredTemperature  <<
//              "     cpu:"<< communic.rank() <<  "         !!!!!!!!!!!!!" << endl ;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LangevinBarostat::registerPython() {
      using namespace espresso::python;

      class_<LangevinBarostat, shared_ptr<LangevinBarostat> >

        ("integrator_LangevinBarostat", init< shared_ptr<System>, shared_ptr<esutil::RNG> >())

        .add_property("gammaP", &LangevinBarostat::getGammaP, &LangevinBarostat::setGammaP)
        .add_property("pressure", &LangevinBarostat::getPressure, &LangevinBarostat::setPressure)
        .add_property("mass", &LangevinBarostat::getMass, &LangevinBarostat::setMass);
    }

  }
}
