#ifndef _INTEGRATOR_HPP
#define _INTEGRATOR_HPP

namespace espresso {
  namespace integrator {

    class integrator {
    protected:
      size_t n_steps;
    public:
      integrator() {};
      integrator(size_t _n_steps): n_steps(_n_steps) {}
      virtual void set_n_steps(size_t _n_steps) {n_steps=_n_steps;}
      virtual size_t get_n_steps(void) const {return n_steps;}
    };

  }
}

#endif
