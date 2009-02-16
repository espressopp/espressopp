#ifndef _ESUTIL_VIRTUAL_FUNCTIONAL_HPP
#define _ESUTIL_VIRTUAL_FUNCTIONAL_HPP

#include <functional>

namespace esutil {
    template <class Argument, class Result>
    class VirtualUnaryFunction
	: public std::unary_function<Argument, Result> {
    public:
	/** General function that is applied to a particle.
	 */
	virtual Result operator()(Argument) = 0;
	virtual ~VirtualUnaryFunction() {}
    };

    template <class Argument1, class Argument2, class Result>
    class VirtualBinaryFunction
	: public std::binary_function<Argument1, Argument2, Result> {
    public:
	/** General function that is applied to a particle.
	 */
	virtual Result operator()(Argument1, Argument2) = 0;
	virtual ~VirtualBinaryFunction() {}
    };
}

#endif
