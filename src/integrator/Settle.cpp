/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "python.hpp"

#include "Settle.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "System.hpp"
#include "bc/BC.hpp"


namespace espressopp {
  using namespace iterator;
  namespace integrator {

    LOG4ESPP_LOGGER(Settle::theLogger, "Settle");

    Settle::Settle(shared_ptr<System> _system, shared_ptr<FixedTupleListAdress> _fixedTupleList,
      real _mO, real _mH, real _distHH, real _distOH)
    : Extension(_system), fixedTupleList(_fixedTupleList),
      mO(_mO), mH(_mH), distHH(_distHH), distOH(_distOH){

        LOG4ESPP_INFO(theLogger, "construct Settle");

        /*
        con1 = integrator->saveOldPos.connect
          (boost::bind(&Settle::saveOldPos, this));
        con2 = integrator->applyConstraints.connect
          (boost::bind(&Settle::applyConstraints, this));
        */

        // initialize settlep variables
        real rmT = 1.0 / (mO + mH + mH);
	mOrmT = mO * rmT;
	mHrmT = mH * rmT;
	real t1 = 0.5 * mO / mH;

	rc = 0.5 * distHH;
	ra = sqrt(distOH * distOH - rc*rc) / (1.0 + t1);
	rb = t1 * ra;
	rra = 1.0 / ra;

        // initialise settlev variables

        mOmH = mO + mH;
        mOmH2 = mOmH * mOmH;
        twicemO = 2*mO;
        twicemH = 2*mH;
        mH2 = mH*mH;

    }

    Settle::~Settle() {
        LOG4ESPP_INFO(theLogger, "~Settle");
        /*
        con1.disconnect();
        con2.disconnect();
        */
    }

    void Settle::disconnect(){
      _befIntP.disconnect();
      _aftIntP.disconnect();
      _aftIntV.disconnect();  // OUT AGAIN?
      _aftIntSlow.disconnect();
    }

    void Settle::connect(){
      // connection to initialisation
      _befIntP  = integrator->befIntP.connect( boost::bind(&Settle::saveOldPos, this));
      _aftIntP  = integrator->aftIntP.connect( boost::bind(&Settle::applyConstraints, this));
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&Settle::correctVelocities, this));   // OUT AGAIN?
      _aftIntSlow  = integrator->aftIntSlow.connect( boost::bind(&Settle::correctVelocities, this));
    }

    void Settle::saveOldPos() {
        oldPos.clear();
        System& system = getSystemRef();
    	// loop over all local molecules
        CellList realCells = system.storage->getRealCells();
        for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            // check if molecule is HHO
            if (molIDs.count(cit->id()) > 0) {

                // lookup cit in tuples, and save AT positions
                FixedTupleListAdress::iterator it;
                it = fixedTupleList->find(&(*cit));

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
    	//System& system = getSystemRef();
        OldPos::iterator it;
        for (it=oldPos.begin(); it != oldPos.end(); it++ ) {
            //std::cout << "in applyConstraints " << system.comm->rank() << " molID " << it->first << std::endl;
            settlep(it->first);
        }
    }

    void Settle::correctVelocities() {

        // call settlev() for every water molecule on node

        System& system = getSystemRef();
        // loop over all local molecules
        CellList realCells = system.storage->getRealCells();
        for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            // check if molecule is HHO
            if (molIDs.count(cit->id()) > 0) {
              settlev(cit->id());
            }
        }
    }


    /*
     * Adapted from NAMD settle1()
     *
     * "This software includes code developed by the Theoretical and Computational
     *  Biophysics Group in the Beckman Institute for Advanced Science and
     *  Technology at the University of Illinois at Urbana-Champaign."
     *
     * Reference for the SETTLE algorithm S. Miyamoto et al.,
     * J. Comp. Chem., 13, 952 (1992).
     *
     */
    void Settle::settlep(longint molID){

    	System& system = getSystemRef();
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
        real dt = integrator->getTimeStep();
        real invdt = 1.0/dt;

    	//std::cout << "\nsettling " << molID << "\n";

    	// positions in previous time step
    	OldPos::const_iterator pos = oldPos.find(molID);
    	if (pos == oldPos.end()) {
    		std::cout << "WARNING: oldPos not found! Skipping settlep().\n";
    		return;
    	}

    	// retrieve particles from tuples
    	Particle* vp = system.storage->lookupRealParticle(molID);
    	FixedTupleListAdress::iterator it;
    	it = fixedTupleList->find(&(*vp));

    	Particle* O  = it->second.at(0);
    	Particle* H1 = it->second.at(1);
    	Particle* H2 = it->second.at(2);

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

        //get unconstrained velocities at v(t+dt)
        Real3D displ1,displ2,displ3;
        bc.getMinimumImageVectorBox(displ1,O->getPos(),pos->second.first); // pos after settle - pos at prev timestep
        Real3D vO=displ1*invdt;
        bc.getMinimumImageVectorBox(displ2,H1->getPos(),pos->second.second);
        Real3D vH1=displ2*invdt;
        bc.getMinimumImageVectorBox(displ3,H2->getPos(),pos->second.third);
        Real3D vH2=displ3*invdt;
        O->setV(vO);
        H1->setV(vH1);
        H2->setV(vH2);

    }

    void Settle::settlev(longint molID){

        //settlev never called, not necessarily debugged

        real dt = integrator->getTimeStep();
        real invdt = 1.0/dt;

        System& system = getSystemRef();
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

        // retrieve particles from tuples
        Particle* vp = system.storage->lookupRealParticle(molID);
        FixedTupleListAdress::iterator it;
        it = fixedTupleList->find(&(*vp));

        Particle* O  = it->second.at(0);
        Particle* H1 = it->second.at(1);
        Particle* H2 = it->second.at(2);

        Real3D vO = O->getV();
        Real3D vH1 = H1->getV();
        Real3D vH2 = H2->getV();

        //get unit vectors along bonds
        Real3D rab;
        Real3D rbc;
        Real3D rca;

        bc.getMinimumImageVectorBox(rab,H1->getPos(),O->getPos());
        bc.getMinimumImageVectorBox(rbc,H2->getPos(),H1->getPos());
        bc.getMinimumImageVectorBox(rca,O->getPos(),H2->getPos());

        real rab2 = rab.sqr();
        real rbc2 = rbc.sqr();
        real rca2 = rca.sqr();

        real rab_abs = sqrt(rab2);
        real rbc_abs = sqrt(rbc2);
        real rca_abs = sqrt(rca2);

        Real3D eab = rab/rab_abs;
        Real3D ebc = rbc/rbc_abs;
        Real3D eca = rca/rca_abs;

        //get relative velocities
        //Real3D vab0r = H1->getV() - O->getV();
        //Real3D vbc0r = H2->getV() - H1->getV();
        //Real3D vca0r = O->getV() - H2->getV();
        Real3D vab0r = vH1 - vO;
        Real3D vbc0r = vH2 - vH1;
        Real3D vca0r = vO - vH2;

        //get components of relative velocities along bonds
        real vab0 = eab * vab0r;
        real vbc0 = ebc * vbc0r;
        real vca0 = eca * vca0r;

        real cosA = (rca2 + rab2 - rbc2)/(2*rca_abs*rab_abs);  //angles depend on constrained geom, no need to recalc, rewrite code to just do it once
        real cosB = (rbc2 + rab2 - rca2)/(2*rbc_abs*rab_abs);  //A=109.527, B=C=35.2819
        real cosC = (rbc2 + rca2 - rab2)/(2*rbc_abs*rca_abs);
        real interm1 = 2*mOmH2 + twicemO*mH*cosA*cosB*cosC - 2*mH2*cosA*cosA - mO*mOmH*(cosB*cosB+cosC*cosC);
        real d = dt * interm1 / (twicemH);

        real interm2 = vab0 * (2*mOmH - mO*cosC*cosC) +
                       vbc0 * (mH*cosC*cosA - mOmH*cosB) +
                       vca0 * (mO*cosB*cosC - twicemH*cosA);
        real tauab = mO*interm2/d;

        real interm3 = vbc0 * (mOmH2 - mH2*cosA*cosA) +
                       vca0 * mO * (mH*cosA*cosB - mOmH*cosC) +
                       vab0 * mO * (mH*cosC*cosA - mOmH*cosB);
        real taubc = interm3/d;

        real interm4 = vca0 * (2*mOmH - mO*cosB*cosB) +
                       vab0 * (mO*cosB*cosC - twicemH*cosA) +
                       vbc0 * (mH*cosA*cosB - mOmH*cosC);
        real tauca = mO*interm4/d;

        Real3D newVa = vO + dt/twicemO*(tauab*eab - tauca*eca);
        Real3D newVb = vH1 + dt/twicemH*(taubc*ebc - tauab*eab);
        Real3D newVc = vH2 + dt/twicemH*(tauca*eca - taubc*ebc);

        O->setV(newVa);
        H1->setV(newVb);
        H2->setV(newVc);

    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Settle::registerPython() {

      using namespace espressopp::python;

      void (Settle::*pyAdd)(longint pid) = &Settle::add;

      class_<Settle, shared_ptr<Settle>, bases<Extension> >
        ("integrator_Settle", init<shared_ptr<System>, shared_ptr<FixedTupleListAdress>, real, real, real, real>())
        .def("add", pyAdd)
        ;
    }
  }
}
