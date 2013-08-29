// ESPP_CLASS
#ifndef _FIXEDSINGLELIST_HPP
#define _FIXEDSINGLELIST_HPP

#include "log4espp.hpp"
#include "python.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <set>
#include <boost/signals2.hpp>


namespace espresso {
	class FixedSingleList : public SingleList{
	  public:
	    typedef std::set<longint> GlobalSingles;

	  protected:
		boost::signals2::connection con1, con2, con3;
		shared_ptr <storage::Storage> storage;
		GlobalSingles globalSingles;
		using SingleList::add;

	  public:
		FixedSingleList(shared_ptr <storage::Storage> _storage);
		virtual ~FixedSingleList();

		/** Add the given particle to the list on this processor if the
		particle belongs to this processor.  Note that
		this routine does not check whether the particle is inserted on
		another processor as well.
		\return whether the particle was inserted on this processor.
		*/
		virtual bool add(longint pid1);
		virtual void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
		void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
		virtual void onParticlesChanged();

	    python::list getSingles();
	    GlobalSingles* getGlobalSingles() {return &globalSingles;};

	    /** Get the number of particles in the GlobalSingle list */
	    int size() {
	    	return globalSingles.size();
	    }

	    static void registerPython();

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

