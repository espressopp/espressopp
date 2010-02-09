#include <iostream>
#include <sys/times.h>
#include <cmath>
#include <cstdlib>

using namespace std;

const int NUM_TESTS = 200000000;

struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

class DirectClass {
public:
  double get(const double x) const {
    return x*x;
  }
};


class VirtualClass {
public:
  virtual double get(const double x) const {
    return sqrt(x);
  }
};

class VirtualSubClass : public VirtualClass {
public:
  virtual double get(const double x) const {
    return x*x;
  }
};


int main() {
  double sum;

  {
    VirtualClass *c = new VirtualSubClass;
    sum = 0.0;
    times(&start_t);
    for (long i = 0; i < NUM_TESTS; i++)
      sum += c->get(drand48());
    times(&stop_t);
    cout << "Did " << NUM_TESTS << " virtual function calls." << endl;
    cout << "user time = " 
	 << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
	 << "s\tsystem time = " 
	 << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
	 << "s" << endl;
    cout << "sum=" << sum << endl;
  }

  {
    DirectClass c;
    sum = 0.0;
    times(&start_t);
    for (long i = 0; i < NUM_TESTS; i++)
      sum += c.get(drand48());
    times(&stop_t);
    cout << "Did " << NUM_TESTS << " direct function calls." << endl;
    cout << "user time = " 
	 << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
	 << "s\tsystem time = " 
	 << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
	 << "s" << endl;
    cout << "sum=" << sum << endl;
  }

  
}
