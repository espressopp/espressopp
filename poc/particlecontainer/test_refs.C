#include <iostream>
#include <vector>
#include <sys/times.h>

using namespace std;

struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

class ThinReference {
    template <class ThinReference, class V> friend class Evaluate;
    int index;
public:
    ThinReference(int i): index(i) {}
};

template <class V>
class FatReference: public ThinReference {
    template <class FatReference, class TV> friend class Evaluate;
    V &v;
public:
    FatReference(int i, V &_v): ThinReference(i), v(_v) {}
};

template<class Reference, class V>
class Evaluate {
    V &v;
public:
    Evaluate(V &_v): v(_v) {}
    virtual void doSomething(Reference ref) {
        v[ref.index] = v[ref.index]*v[ref.index];
    }
};

int main()
{
    const int NUM_TESTS = 2000;
    const int VSIZE = 100000;

    vector<int> vec(VSIZE, 22);

    // benchmark thin
    for (int repeat = 0; repeat < 5; ++repeat) {
        Evaluate<ThinReference, vector<int> > *eval = new Evaluate<ThinReference, vector<int> >(vec);
        std::cout << "Testing thin..." << std::endl;
        times(&start_t);
        {
            for (int r=0; r < NUM_TESTS; r++) {
                for (int i=0; i < VSIZE; i++) {
                    eval->doSomething(ThinReference(i));
                }
            }
        }
        times(&stop_t);
        std::cout << "user time = " 
                  << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
                  << "s\tsystem time = "
                  << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
                  << "s" << std::endl;
    }

    // benchmark fat
    for (int repeat = 0; repeat < 5; ++repeat) {
        Evaluate<FatReference<vector<int> >, vector<int> > *eval =
            new Evaluate<FatReference<vector<int> >, vector<int> >(vec);
        std::cout << "Testing fat..." << std::endl;
        times(&start_t);
        {
            for (int r=0; r < NUM_TESTS; r++) {
                for (int i=0; i < VSIZE; i++) {
                    eval->doSomething(FatReference<vector<int> >(i, vec));
                }
            }
        }
        times(&stop_t);
        std::cout << "user time = " 
                  << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
                  << "s\tsystem time = "
                  << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
                  << "s" << std::endl;
    }

    // benchmark thin reference
    for (int repeat = 0; repeat < 5; ++repeat) {
        Evaluate<const ThinReference &, vector<int> > *eval =
            new Evaluate<const ThinReference &, vector<int> >(vec);
        std::cout << "Testing thin referenced..." << std::endl;
        times(&start_t);
        {
            for (int r=0; r < NUM_TESTS; r++) {
                for (int i=0; i < VSIZE; i++) {
                    eval->doSomething(ThinReference(i));
                }
            }
        }
        times(&stop_t);
        std::cout << "user time = " 
                  << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
                  << "s\tsystem time = "
                  << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
                  << "s" << std::endl;
    }

    // benchmark fat
    for (int repeat = 0; repeat < 5; ++repeat) {
        Evaluate<const FatReference<vector<int> > &, vector<int> > *eval =
            new Evaluate<const FatReference<vector<int> > &, vector<int> >(vec);
        std::cout << "Testing fat referenced..." << std::endl;
        times(&start_t);
        {
            for (int r=0; r < NUM_TESTS; r++) {
                for (int i=0; i < VSIZE; i++) {
                    eval->doSomething(FatReference<vector<int> >(i, vec));
                }
            }
        }
        times(&stop_t);
        std::cout << "user time = " 
                  << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
                  << "s\tsystem time = "
                  << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
                  << "s" << std::endl;
    }

    return 0;
}

/*
  Local Variables:
  compile-command: "g++ -Wall -static -O3 test_refs.C -o test_refs && ./test_refs"
  End:
*/
