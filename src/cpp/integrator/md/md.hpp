#ifndef _MD_INTEGRATOR_HPP
#define _MD_INTEGRATOR_HPP

#include "integrator/integrator.hpp"

namespace espresso {
  namespace integrator {

    class md: public integrator {
    protected:
      real md_time;
      real time_step;
      real time_step_Sqr;
    public:
      md() {};
      md(real _time_step, real _md_time): time_step(_time_step), md_time(_md_time) {};
      virtual void set_time_step(real _time_step) {
        time_step = _time_step;
        time_step_Sqr = time_step * time_step;
      }
      virtual void set_md_time(real _md_time) {md_time=_md_time;}
      virtual real get_time_step(void) const {return time_step;}
      virtual real get_time_step_Sqr(void) const {return time_step_Sqr;}
      virtual real get_md_time(void) const {return md_time;}
    };

  }
}

#endif
