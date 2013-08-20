#include "python.hpp"
#include "OrderParameter.hpp"

namespace espresso {
  namespace analysis {


    dcomplex OrderParameter::SphHarm(int l, int m, Real3D r) {
      real d = r.abs(); // distance between two particles
      real theta, phi;

      // defining theta and phi
      theta = acos( r[2] / d );   // in radians
      
      // phi is defined as http://en.wikipedia.org/wiki/Atan2
      // problem is x = y = 0, we will define it like 0
      if(r[0] > 0.0){
        phi = atan( r[1]/r[0] );
      }
      else if( r[0] < 0.0 && r[1] >= 0.0 ){
        phi = atan( r[1]/r[0] ) + M_PIl;
      }
      else if( r[0] < 0.0 && r[1] < 0.0 ){
        phi = atan( r[1]/r[0] ) - M_PIl;
      }
      else if( r[0] == 0.0 && r[1] > 0.0 ){
        phi = M_PIl;
      }
      else if( r[0] == 0.0 && r[1] < 0.0 ){
        phi = -M_PIl;
      }
      else{
        // x = 0; y = 0;
        phi = 0.0;
      }
      
      dcomplex res = boost::math::spherical_harmonic(l, m, theta, phi);
      
      //std::cout << " boost gets: l=" << l << "  m="<<m<<"  theta="<<theta<<"  phi="<<phi << std::endl;
      //std::cout << " SphHarm returns: " << res << std::endl;
      return res;
    }
    
    void OrderParameter::registerPython() {
      using namespace espresso::python;
      class_<OrderParameter, bases< AnalysisBase > >
        ("analysis_OrderParameter", init< shared_ptr<System>, real, int, real >())
        .add_property("l", 
              &OrderParameter::getAngularMomentum,
              &OrderParameter::setAngularMomentum)
        .add_property("cutoff", 
              &OrderParameter::getCutoff,
              &OrderParameter::setCutoff)
        .add_property("threshold", 
              &OrderParameter::getThreshold,
              &OrderParameter::setThreshold)
      ;
    }
  }
}
