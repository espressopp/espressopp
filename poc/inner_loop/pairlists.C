#include <iostream>
#include <sys/times.h>
#include <vector>
#include <utility>
#include <stdlib.h>

using namespace std;

const int N = 10000;
const int N_PARTNERS = 100;
const int NUM_TESTS = 1000;

struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

typedef pair<double*,double*> Pair;

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
  
  // for the pair implementation
  vector< Pair > pairs;
  vector< Pair >::const_iterator pit;
  
  // for the loop implementation
  vector< vector< double* > > vl;
  vl.resize(N);
  
  for (i = 0; i< N; i++) {
    p1 = &v[i];
    for (j = 0; j < N_PARTNERS; j++) {
      r = static_cast<int>(drand48()*N);
      p2 = &v[r];
      // fill the verlet list
      vl[i].push_back(p2);
      // push a new pair
      pairs.push_back(make_pair(p1, p2));
    }
  }
  
  cout << "Created " << pairs.size() << " pairs." << endl;
  
  // ****************************************
  // benchmark pairs
  // ****************************************
  sum = 0.0;
  count = 0;
  times(&start_t);
  for (i=0; i< NUM_TESTS; i++)
    for (pit = pairs.begin(); pit != pairs.end(); pit++) {
      sum += *(pit->first) - *(pit->second);
      count++;
    }
  times(&stop_t);
  cout << "Pairs implementations:" << endl;
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = " 
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;
  cout << "sum=" << sum << "\tcount=" << count << endl;

  // ****************************************
  // benchmark Verlet list
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
	count++;
      }
    }
  times(&stop_t);
  cout << "VL implementations:" << endl;
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = " 
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;
  cout << "sum=" << sum << "\tcount=" << count << endl;
  
  

}
