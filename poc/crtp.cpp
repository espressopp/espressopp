#include <iostream>

/*
  This is a test program to test the so called
  curiously recurrent template pattern (CRTP), and
  to demonstrate the powers of the static_cast<>
  in comparison to a simple C-style cast.

  The class Derived is the central class; it consists of a
  Spacer class, the class following the CRTP and some data
  of its own. The Spacer just occupies some space so that
  the base address of the CRTP object differs from the
  base address of its Derived parent. As a consequence, the
  static_cast from CRTP to Derived actually needs to _subtract_
  something from the memory address.  Obtaining similar behaviour
  in C would require some serious pointer arithmetics.
*/

// calculate bytewise relative addresses for comparison
ptrdiff_t relAddr(void *p, void *b) {
    return static_cast<char *>(p) - static_cast<char *>(b);
}

struct Spacer {
    int spacerVal;
};

template<class Derived>
class CRTP {
public:
    int crtpVal;

    int getValueIndirect() {
        Derived *derived = static_cast<Derived *>(this);
        std::cout << "CRTP: my Derived at: " <<  relAddr(derived, this) << std::endl;
        return derived->getValue();
    }
};

class Derived: public Spacer, public CRTP<Derived> {
    int derivedVal;    

    ptrdiff_t relAddr(void *p) { return ::relAddr(p, this); }

public:
    Derived(): derivedVal(42) {
        std::cout << "my spacer val at:  " <<  relAddr(&spacerVal) << std::endl
                  << "my crtp val at:    " <<  relAddr(&crtpVal) << std::endl
                  << "my derived val at: " <<  relAddr(&derivedVal) << std::endl
                  << "my spacer at: " << relAddr(static_cast<Spacer *>(this)) << std::endl
                  << "my crtp at:   " << relAddr(static_cast<CRTP<Derived> *>(this)) << std::endl;
    }

    int getValue() { return derivedVal; }
};

int main()
{
    Derived test;
    std::cout << "-----> direct:" << test.getValue() << " indirect: " << test.getValueIndirect() << std::endl;
}
