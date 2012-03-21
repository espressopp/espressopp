// ESPP_CLASS
#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"
#include "mpi.hpp"
#include "python.hpp"
#include "boost/enable_shared_from_this.hpp"
#include "interaction/Interaction.hpp"

namespace espresso {

  namespace esutil { class RNG; }

  class System : public enable_shared_from_this< System > {

  public:

    System();
    System(python::object _pyobj);

    shared_ptr< mpi::communicator > comm;

    shared_ptr< storage::Storage > storage;
    shared_ptr< bc::BC > bc;
    shared_ptr< esutil::RNG > rng;

    interaction::InteractionList shortRangeInteractions;

    real skin;  //<! skin used for VerletList
    
    real maxCutoff;     // maximal cutoff over all of the interactions

    bool CommunicatorIsInitialized;

    shared_ptr< System > getShared() { 
      return shared_from_this();
    }
    
    void scaleVolume(real s, bool particleCoordinates);
    void scaleVolume(Real3D s, bool particleCoordinates);
    void scaleVolume3D(Real3D s);
    void setTrace(bool flag);
    void addInteraction(shared_ptr< interaction::Interaction > ia);
    void removeInteraction(int i);
    shared_ptr< interaction::Interaction > getInteraction(int i);
    int getNumberOfInteractions();
    static void registerPython();

  };
}
#endif
