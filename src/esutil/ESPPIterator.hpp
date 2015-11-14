/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
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

// ESPP_CLASS
#ifndef _ESUTIL_ESPPITERATOR_HPP
#define _ESUTIL_ESPPITERATOR_HPP

namespace espressopp {
  namespace esutil {
    template < class STLContainer > 
    class ESPPIterator {
      typedef typename STLContainer::value_type value_type;
    public:
      ESPPIterator()
	: stlIt(), stlEnd()
      {}

      ESPPIterator(STLContainer &container) 
	: stlIt(container.begin()), stlEnd(container.end())
      {}
      
      ESPPIterator(typename STLContainer::iterator begin, 
		   typename STLContainer::iterator end) 
	: stlIt(begin), stlEnd(end)
      {}
      
      ESPPIterator &operator++() { ++stlIt; return *this; }
      bool isValid() const { return stlIt != stlEnd; }
      bool isDone() const { return !isValid(); }
      
      value_type &operator*() const { return *stlIt; }
      value_type *operator->() const { return &**this; }

      typename STLContainer::iterator 
      getSTLIterator() { return stlIt; }
      
    private:
      typename STLContainer::iterator stlIt;
      typename STLContainer::iterator stlEnd;
    };

    template < class STLContainer >
    class ESPPContainer : public STLContainer {
    protected:
      typedef STLContainer Super;
    public:
      typedef ESPPIterator<STLContainer> Iterator;
    };
  }
}

#endif
