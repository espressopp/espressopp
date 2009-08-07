#ifndef _ESUTIL_MULTISIGNALCONNECTIONS_HPP
#define _ESUTIL_MULTISIGNALCONNECTIONS_HPP

#include "types.hpp"
#include <boost/signals2.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <map>
#include <utility>

namespace espresso {
  namespace esutil {
    // manages multiple connections
    class MultiSignalConnections 
      : public std::multimap< void*, boost::signals2::connection >  {

    public:
      template < typename SignalType, typename SharedPtr, typename MethodPtr >
      void add(SignalType &signal,
	       SharedPtr obj, MethodPtr method, 
	       boost::signals2::connect_position at 
	       = boost::signals2::at_back) {
	boost::signals2::connection connection 
	  = signal.connect(at, boost::bind(method, obj, _1));
	insert(std::make_pair(static_cast< void* >(obj.get()), connection));
      }

      // disconnect and remove the connections of this listener
      template < typename SharedPtr >
      void remove(SharedPtr listener) {
	std::pair< iterator, iterator >
	  range = equal_range(static_cast< void* >(listener.get())); 
	BOOST_FOREACH(value_type connection, range)
	  {
	    connection.second.disconnect();
	  }
	erase(range.first, range.second);
      }

      ~MultiSignalConnections() {
	BOOST_FOREACH(value_type connection, *this) {
	  connection.second.disconnect();
	}
	erase(begin(), end());
      }
    };
  }
} 

#endif
