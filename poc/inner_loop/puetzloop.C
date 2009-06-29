#include <iostream>
#include <vector>
#include <sys/times.h>
#include <cmath>

using namespace std;

const int N = 10000;
const int N_PARTNERS = 100;
const int NUM_TESTS = 1000;

struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

int main() {
  int i, j, k;
  int count;
  double sum;
  vector< double > v;
  vector< double >::iterator vit;
  for (i=0; i< N; i++)
    v.push_back(drand48());
  
  cout << "Created " << v.size() << " random doubles." << endl;

  double *p1;
  double *p2;
  int r;

  // for the loop implementation
  vector< vector< double* > > vl;
  vl.resize(N);
  
  for (i = 0; i< N; i++) {
    p1 = &v[i];
    for (j = 0; j < N_PARTNERS; j++) {
      r = static_cast<int>(drand48() * N);
      p2 = &v[r];
      // fill the verlet list
      vl[i].push_back(p2);
    }
  }
  
  cout << "Created " << (N*N_PARTNERS) << " interactions." << endl;

  {
    // ****************************************
    // Benchmark Verlet list
    // ****************************************
    vector< double* >::const_iterator it, endit;
    register double x1;
    sum = 0.0;
    count = 0;
    times(&start_t);
    for (i = 0; i < NUM_TESTS; i++)
      for (j=0; j < N; j++) {
	x1 = v[j];
	endit = vl[j].end();
	for (it = vl[j].begin(); it != endit; it++) {
	  sum += x1 - **it;
	}
      }
    times(&stop_t);
    cout << "VL implementations:" << endl;
    cout << "user time = " 
	 << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
	 << "s\tsystem time = " 
	 << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
	 << "s" << endl;
    cout << "sum=" << sum << endl;
  }

  {
    // ****************************************
    // Benchmark Puetzified Verlet list
    // ****************************************
    vector< double* >::const_iterator it, endit;
    register double x1;
    sum = 0.0;
    times(&start_t);
    for (i = 0; i < NUM_TESTS; i++)
      for (j=0; j < N; j++) {
	x1 = v[j];
	endit = vl[j].end();
	for (it = vl[j].begin(); it != endit; it++) {
	  sum += x1 - **it;
	}
      }
    times(&stop_t);
    cout << "VL implementations:" << endl;
    cout << "user time = " 
	 << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
	 << "s\tsystem time = " 
	 << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
	 << "s" << endl;
    cout << "sum=" << sum << endl;
  }

}
