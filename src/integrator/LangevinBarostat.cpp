#include "python.hpp"
#include "LangevinBarostat.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

#include "interaction/Interaction.hpp"
#include "esutil/RNG.hpp"

#include "bc/BC.hpp"

#include "mpi.hpp"

namespace espresso {

  using namespace std;
  
  using namespace iterator;
  using namespace interaction;

  namespace integrator {

    LOG4ESPP_LOGGER(LangevinBarostat::theLogger, "LangevinBarostat");

    // TODO It is very sensitive to parameters (the ensemble). The manual how to choose the
    // parameters for simple system should be written.
    LangevinBarostat::LangevinBarostat(shared_ptr<System> _system, shared_ptr< esutil::RNG > _rng) : SystemAccess(_system), rng(_rng){
      // external parameters
      gammaP = 0.0;
      mass = 0.0;
      externalPressure = 0.0;
      
      // local variable
      momentum = 0.0;
      momentum_mass = 0.0;

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

    void LangevinBarostat::updVolume(real dt){
      System& system = getSystemRef();
      
      // The volume is scaled according to the equations V(t+1/2*dt) = V(t) + 1/2*dt*V'; V' = d*V*pe/W
      // in order to accelerate  1.5 * dt := 0.5dt * 3
      real volScale = 1 + dt * 1.5 * momentum_mass;
      real scale_factor = 0;
      if( volScale < 0.0){
        cout << "The volume scaling is not possible. Volume scaling factor: "<< volScale << endl;
        // TODO error handling
        exit(-1);
      }
      else
        scale_factor = pow( volScale, 1./3.);  // calculating the current scaling parameter
      
      system.scaleVolume( scale_factor, false);
    }
    
    // TODO should be optimized!!!
    /*
     *  TODO now it is valid for nonbonded systems because the value of degreees of freedom
     *  Nf = 3*N, N - number of particles, 3 - d-dimensional system (d=3). Thus d/Nf is 
     *  replaced by 1/N.
     */
    void LangevinBarostat::updVolumeMomentum(real dt_2){
      // dt_2 is timestep/2. 
      System& system = getSystemRef();
      Real3D Li = system.bc -> getBoxL(); // getting the system size
      real V = Li[0] * Li[1] * Li[2];     // system volume
      real m3V = 3 * V;
      
      // get a random value and distribute the same value over all of the CPUs
      mpi::communicator communic = *system.comm;
      real rannum; 
      if (communic.rank() == 0) rannum = (*rng)()-0.5; //rannum = rng->normal();
      mpi::broadcast(communic, rannum, 0);
      
      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      // it's not efficient to use the analysis.Pressure because of double calculation of m*v*v
      real v2sum;
      real v2 = 0.0;
      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        const Particle& p = *cit;
        v2 += p.mass() * (p.velocity() * p.velocity());
      }
      mpi::all_reduce( communic, v2, v2sum, std::plus<real>());
      real p_kinetic = v2sum;
      
      // compute the short-range nonbonded contribution
      real rij_dot_Fij = 0.0;
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        rij_dot_Fij += srIL[j]->computeVirial();
      }
      real p_nonbonded = 0.0;
      mpi::all_reduce( communic, rij_dot_Fij, p_nonbonded, std::plus<real>());
      
      // @TODO one should check that X is not the instantaneous pressure, 
      // because it does not include the white noise from thermostat
      real X = p_kinetic + p_nonbonded;
      /*
       * in order not to do double calculations X is not divided by 3V and term in the next line
       * 3V * (X - externalPressure) is replaced by (X - 3V * externalPressure)
       * local momentum derivative; in order to optimize term pref6 * v2sum could be coupled with X
       */
      real pe_deriv = (X - m3V * externalPressure) + pref6 * v2sum + pref4 * momentum + pref5 * rannum;
      
      momentum += dt_2 * pe_deriv;
      momentum_mass = momentum / mass; // momentum is normalized by mass already
    }
    
    void LangevinBarostat::updForces(){
      LOG4ESPP_DEBUG(theLogger, "barostating");

      System& system = getSystemRef();
      CellList cells = system.storage->getRealCells();

      real factor = pref3  * momentum_mass;
      for(CellListIterator cit(cells); !cit.isDone(); ++cit){
        frictionBarostat(*cit, factor);
      }
    }
    
    real LangevinBarostat::updDisplacement(){
      // should return pe/W
      return (momentum_mass); // momentum is normalized by fictitious mass already
    }

    void LangevinBarostat::frictionBarostat(Particle& p, real factor){
      p.force() += factor * p.velocity() * p.mass();
      
      LOG4ESPP_TRACE(theLogger, "new particle force = " << p.force());
    }

    void LangevinBarostat::initialize(real timestep, real desiredTemperature){
      // calculate the prefactors
      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
                               ", gammaP = " << gammaP << 
                               ", external pressure = " << externalPressure <<
                               ", fictitious mass = " << mass );

      // TODO implement Nf for the system with constarins
      /* Nf - degrees of freedom. For N particles in d-dimentional system without constrains
       * Nf = d * N
       */
      // determine the number of local particles and total particles
      System& system = getSystemRef();
      int Nsum = 0;
      int N = system.storage->getNRealParticles();
      boost::mpi::all_reduce(*mpiWorld, N, Nsum, std::plus<int>());
      
      pref6 = 1./(double)Nsum;
      
      pref3 = - (1+pref6);
      
      // pressure friction prefactor
      pref4 = -gammaP;
      
      // uniform distribution prefactor. (it could be used instead of normal distribution)
      pref5 = 2.0 * sqrt(2.0 * desiredTemperature * gammaP * mass / timestep);
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
