#ifndef _ESUTIL_SIGNALBINDINGTRAITS_HPP
#define _ESUTIL_SIGNALBINDINGTRAITS_HPP

#include "types.hpp"
#include <boost/signals2.hpp>
#include <boost/bind.hpp>

namespace espresso {
  namespace esutil {
    /// automatic binding for boost::signals2
    template < typename SignalType >
    struct SignalBindingTraits {};

#define SIGNALBINDINGTRAITS_REPEATER(z, i, what) , what##i
#define SIGNALBINDINGTRAITS_REPEAT(cnt, what)                           \
    BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_INC(cnt), SIGNALBINDINGTRAITS_REPEATER, what)
#define SIGNALBINDINGTRAITS_CLASS(cnt)                                  \
    boost::signals2::signal ## cnt<Result                               \
    SIGNALBINDINGTRAITS_REPEAT(cnt, Arg) >                              \

#define SIGNALBINDINGTRAITS(z, cnt, data)                               \
    template < typename Result                                          \
               SIGNALBINDINGTRAITS_REPEAT(cnt, typename Arg) >          \
    struct SignalBindingTraits< typename                                \
                                SIGNALBINDINGTRAITS_CLASS(cnt) >        \
    {                                                                   \
      template< typename MethodPtr, typename ObjectPtr > static         \
        typename SIGNALBINDINGTRAITS_CLASS(cnt)::slot_type              \
        bind(MethodPtr method, ObjectPtr obj) {                         \
        return boost::bind(method, obj                                  \
                           SIGNALBINDINGTRAITS_REPEAT(cnt, _));         \
      }                                                                 \
    };

    BOOST_PP_REPEAT(BOOST_SIGNALS2_MAX_ARGS, SIGNALBINDINGTRAITS, );
    
#undef SIGNALBINDINGTRAITS_REPEATER
#undef SIGNALBINDINGTRAITS_REPEAT
#undef SIGNALBINDINGTRAITS_CLASS
#undef SIGNALBINDINGTRAITS

  }
}
#endif
