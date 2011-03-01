// ESPP_CLASS
#ifndef _ESUTIL_ESPPITERATOR_HPP
#define _ESUTIL_ESPPITERATOR_HPP

namespace espresso {
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
