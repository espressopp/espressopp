#include "python.hpp"

#include "Settle.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "System.hpp"


namespace espresso {

    LOG4ESPP_LOGGER(Settle::theLogger, "Settle");

    Settle::Settle(shared_ptr<storage::Storage> _storage,
    		shared_ptr<integrator::VelocityVerlet> _integrator,
    		real _mO, real _mH, real _distHH, real _distOH)
    : storage(_storage), integrator(_integrator),
      mO(_mO), mH(_mH), distHH(_distHH), distOH(_distOH){

        LOG4ESPP_INFO(theLogger, "construct Settle");

        /*
        con1 = integrator->saveOldPos.connect
          (boost::bind(&Settle::saveOldPos, this));
        con2 = integrator->applyConstraints.connect
          (boost::bind(&Settle::applyConstraints, this));
        */

        // initialize settle variables
        real rmT = 1.0 / (mO + mH + mH);
		mOrmT = mO * rmT;
		mHrmT = mH * rmT;
		real t1 = 0.5 * mO / mH;

		rc = 0.5 * distHH;
		ra = sqrt(distOH * distOH - rc*rc) / (1.0 + t1);
		rb = t1 * ra;
		rra = 1.0 / ra;

		fixedtupleList = storage->getFixedTuples();

    }

    Settle::~Settle() {
        LOG4ESPP_INFO(theLogger, "~Settle");
        /*
        con1.disconnect();
        con2.disconnect();
        */
    }


    void Settle::saveOldPos() {
        oldPos.clear();

    	// loop over all local molecules
        CellList realCells = storage->getRealCells();
        for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            // check if molecule is HHO
            if (molIDs.count(cit->id()) > 0) {

                // lookup cit in tuples, and save AT positions
                FixedTupleList::iterator it;
                it = fixedtupleList->find(&(*cit));

                oldPos.insert(std::make_pair(cit->id(),
                        Triple<Real3D,Real3D,Real3D>(
                                (it->second.at(0))->getPos(),
                                (it->second.at(1))->getPos(),
                                (it->second.at(2))->getPos())));
            }
        }
    }


    void Settle::applyConstraints() {

        // call settlep() for every water molecule on node
        OldPos::iterator it;
        for (it=oldPos.begin(); it != oldPos.end(); it++ ) {
            // settle it down!
            settlep(it->first);
        }
	}


    /*
     * Adapted from NAMD settle1()
     *
     * Reference for the SETTLE algorithm S. Miyamoto et al.,
     * J. Comp. Chem., 13, 952 (1992).
     *
     */
    void Settle::settlep(longint molID){

    	//std::cout << "\nsettling " << molID << "\n";

    	// positions in previous time step
    	OldPos::const_iterator pos = oldPos.find(molID);
    	if (pos == oldPos.end()) {
    		std::cout << "WARNING: oldPos not found! Skipping settle().\n";
    		return;
    	}


    	// retrieve particles from tuples
    	Particle* vp = storage->lookupRealParticle(molID);
    	FixedTupleList::iterator it;
    	it = fixedtupleList->find(&(*vp));

    	Particle* O  = it->second.at(0);
    	Particle* H1 = it->second.at(1);
    	Particle* H2 = it->second.at(2);

    	/*
    	std::cout << "old  O pos: " << O->getPos() << "\n";
    	std::cout << "old H1 pos: " << H1->getPos() << "\n";
    	std::cout << "old H2 pos: " << H2->getPos() << "\n";
    	*/


    	// --- Step1 A1' ---
    	// vectors in the plane of the original positions
    	// previous positions OHH
    	Real3D b0 = pos->second.second - pos->second.first; // H1.pos - O.pos
    	Real3D c0 = pos->second.third - pos->second.first; // H2.pos - O.pos

    	// new center of mass
    	// present positions OHH
    	Real3D d0 = O->position() * mOrmT + ((H1->position() + H2->position())*mHrmT);

    	Real3D a1 = O->position() - d0;
    	Real3D b1 = H1->position() - d0;
    	Real3D c1 = H2->position() - d0;

    	// Vectors describing transformation from original coordinate system to
    	// the 'primed' coordinate system
    	Real3D n0 = b0.cross(c0);
    	Real3D n1 = a1.cross(n0);
    	Real3D n2 = n0.cross(n1);

    	// unit vectors
    	n0 = n0/n0.abs();
    	n1 = n1/n1.abs();
    	n2 = n2/n2.abs();

    	b0 = Real3D(n1*b0, n2*b0, n0*b0); // note: b0.z is never referenced again
    	c0 = Real3D(n1*c0, n2*c0, n0*c0); // note: c0.z is never referenced again

    	// shuffle the order of operations around from the normal algorithm so
    	// that we can double pump sqrt() with n2 and cosphi at the same time
    	// these components are usually computed down in the canonical water block
    	real A1Z = n0 * a1;
    	b1 = Real3D(n1*b1, n2*b1, n0*b1);
    	c1 = Real3D(n1*c1, n2*c1, n0*c1);


    	// --- Step2 A2' ---
    	// now we can compute positions of canonical water
    	real sinphi = A1Z * rra;
    	real tmp = 1.0 - sinphi*sinphi;
    	real cosphi = sqrt(tmp);
    	real sinpsi = (b1[2] - c1[2]) / (2.0 * rc * cosphi);
    	tmp = 1.0 - sinpsi*sinpsi;
    	real cospsi = sqrt(tmp);

    	real rbphi = -rb * cosphi;
    	real tmp1 = rc * sinpsi*sinphi;
    	real tmp2 = rc * sinpsi*cosphi;

    	Real3D a2(0, ra * cosphi, ra * sinphi);
    	Real3D b2(-rc * cospsi, rbphi - tmp1, -rb * sinphi + tmp2);
    	Real3D c2( rc * cosphi, rbphi + tmp1, -rb * sinphi - tmp2);


    	// --- Step3 al, be, ga ---
    	// there are no a0 terms because we've already subtracted the term off
    	// when we first defined b0 and c0.
    	real alpha = b2[0] * (b0[0] - c0[0]) + b0[1] * b2[1] + c0[1] * c2[1];
    	real beta  = b2[0] * (c0[1] - b0[1]) + b0[0] * b2[1] + c0[0] * c2[1];
    	real gama  = b0[0] * b1[1] - b1[0] * b0[1] + c0[0] * c1[1] - c1[0] * c0[1];

    	real a2b2 = alpha*alpha + beta*beta;
    	real sintheta = (alpha*gama - beta*sqrt(a2b2 - gama*gama))/a2b2;


    	// --- Step4 A3' ---
    	real costheta = sqrt(1.0 - sintheta*sintheta);

    	Real3D a3(-a2[1] * sintheta,
    			a2[1] * costheta,
    			A1Z);

    	Real3D b3(b2[0] * costheta - b2[1] * sintheta,
    			b2[0] * sintheta + b2[1] * costheta,
    			b1[2]);

    	Real3D c3(-b2[0] * costheta - c2[1] * sintheta,
    			-b2[0] * sintheta + c2[1] * costheta,
    			c1[2]);


    	// --- Step5 A3 ---
    	// undo the transformation; generate new normal vectors from the transpose.
    	Real3D m1(n1[0], n2[0], n0[0]);
    	Real3D m2(n1[1], n2[1], n0[1]);
    	Real3D m0(n1[2], n2[2], n0[2]);

    	// new positions
    	O->position() = Real3D(a3*m1, a3*m2, a3*m0) + d0;
    	H1->position() = Real3D(b3*m1, b3*m2, b3*m0) + d0;
    	H2->position() = Real3D(c3*m1, c3*m2, c3*m0) + d0;

    	/*
    	std::cout << " new  O pos: " << O->getPos() << "\n";
        std::cout << " new H1 pos: " << H1->getPos() << "\n";
        std::cout << " new H2 pos: " << H2->getPos() << "\n";
        */
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Settle::registerPython() {

      using namespace espresso::python;

      void (Settle::*pyAdd)(longint pid) = &Settle::add;

      class_<Settle, shared_ptr<Settle> >
        ("Settle", init<
        		shared_ptr<storage::Storage>,
        		shared_ptr<integrator::VelocityVerlet>,
        		real, real, real, real>())
        .def("add", pyAdd)
        ;
    }

}
