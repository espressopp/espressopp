#ifndef _VIRTUAL_FUNCTIONAL_HPP
#define _VIRTUAL_FUNCTIONAL_HPP

#include <functional>

namespace util {
    template <class Result, class Argument>
    class VirtualUnaryFunction
	: public std::unary_function<Result, Argument> {
    public:
	/** General function that is applied to a particle.
	 */
	virtual Result operator()(Argument) = 0;
    };

    template <class Result, class Argument1, class Argument2>
    class VirtualBinaryFunction
	: public std::binary_function<Result, Argument1, Argument2> {
    public:
	/** General function that is applied to a particle.
	 */
	virtual Result operator()(Argument1, Argument2) = 0;
    };
}

#endif
