/*
 Copyright (C) 2012-2016 Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
 
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
#include "LJcos.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
	namespace interaction {
		
		typedef class VerletListInteractionTemplate <LJcos>
		VerletListLJcos;
		typedef class VerletListAdressInteractionTemplate <LJcos, Tabulated>
		VerletListAdressLJcos;
		typedef class VerletListHadressInteractionTemplate <LJcos, Tabulated>
		VerletListHadressLJcos;
		typedef class CellListAllPairsInteractionTemplate <LJcos>
		CellListLJcos;
		typedef class FixedPairListInteractionTemplate <LJcos>
		FixedPairListLJcos;
		
		//////////////////////////////////////////////////
		// REGISTRATION WITH PYTHON
		//////////////////////////////////////////////////
		void LJcos::registerPython() {
			using namespace espressopp::python;
			
			class_< LJcos, bases< Potential > >
			("interaction_LJcos", init< real >())
			.add_property("phi", &LJcos::getPhi, &LJcos::setPhi)
			.add_property("sigma", &LJcos::getSigma, &LJcos::setSigma)
			.def_pickle(LJcos_pickle())
			;
			
			class_< VerletListLJcos, bases< Interaction > >
			("interaction_VerletListLJcos", init< shared_ptr<VerletList> >())
			.def("getVerletList", &VerletListLJcos::getVerletList)
			.def("setPotential", &VerletListLJcos::setPotential, return_value_policy< reference_existing_object >())
			.def("getPotential", &VerletListLJcos::getPotential, return_value_policy< reference_existing_object >())
			;
			
			class_< VerletListAdressLJcos, bases< Interaction > >
			("interaction_VerletListAdressLJcos",
			 init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
			.def("setPotentialAT", &VerletListAdressLJcos::setPotentialAT)
			.def("setPotentialCG", &VerletListAdressLJcos::setPotentialCG);
			;
			
			class_< VerletListHadressLJcos, bases< Interaction > >
			("interaction_VerletListHadressLJcos",
			 init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
			.def("setPotentialAT", &VerletListHadressLJcos::setPotentialAT)
			.def("setPotentialCG", &VerletListHadressLJcos::setPotentialCG);
			;
			
			class_< CellListLJcos, bases< Interaction > >
			("interaction_CellListLJcos", init< shared_ptr< storage::Storage > >())
			.def("setPotential", &CellListLJcos::setPotential);
	  ;
			
			class_< FixedPairListLJcos, bases< Interaction > >
			("interaction_FixedPairListLJcos",
			 init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LJcos> >())
			.def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LJcos> >())
			.def("setPotential", &FixedPairListLJcos::setPotential)
			.def("setFixedPairList", &FixedPairListLJcos::setFixedPairList)
			.def("getFixedPairList", &FixedPairListLJcos::getFixedPairList)
			;
		}
		
	}
}
