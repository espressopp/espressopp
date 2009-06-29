#include <iostream>
#include <vector>
#include <sys/times.h>

/* this code just demonstrates that a class Integer that contains just an
   int and behaves like one is also treated as one by the compiler. This
   means that one can use a specific index class for an array where the
   generator routine assures that the index is always in range */

const int NUM_TESTS = 10000;
const int VSIZE = 100000;

struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

void tagit() __attribute__((noinline));

void tagit() { static int c = 0; c++; }

class Integer {
    int index;

public:
    Integer(int i): index(i) {}
    operator const int () const { return index; }
};

class OtherInteger {
    friend class Integer;
    int index;

public:
    OtherInteger(int i): index(i) {}
    operator Integer() const {
        tagit();
        return Integer(index);
    }
};

template<class index>
class Property {
    std::vector<double> data;

public:
    Property() {
        data.resize(VSIZE);
    }

    double &operator[](index i) {
        return data[i];
    }
};


template<class t, class altt>
double check(Property<t> &p)
    __attribute__((noinline));

template<class t, class altt>
double check(Property<t> &p) {
    double sum = 0;
    for (int i = 0; i < VSIZE; ++i)
        sum += p[altt(i)];
 
    return sum;
}

int main()
{
    // benchmark indexing with an integer
    std::cout << "Testing integer..." << std::endl;
    {
        Property<int> pint;
        times(&start_t);
        for (int i=0; i < NUM_TESTS; i++) {
            check<int,int>(pint);
        }
        times(&stop_t);
    }
    std::cout << "user time = " 
              << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
              << "s\tsystem time = "
              << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
              << "s" << std::endl;
    
    // benchmark indexing with an integer
    std::cout << "Testing integer..." << std::endl;
    {
        Property<int> pint;
        times(&start_t);
        for (int i=0; i < NUM_TESTS; i++) {
            check<int,int>(pint);
        }
        times(&stop_t);
    }
    std::cout << "user time = " 
              << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
              << "s\tsystem time = "
              << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
              << "s" << std::endl;
    
    // benchmark indexing with a wrapped integer
    std::cout << "Testing wrapped integer..." << std::endl;
    {
        Property<Integer> pint;
        times(&start_t);
        for (int i=0; i < NUM_TESTS; i++) {
            check<Integer,Integer>(pint);
        }
        times(&stop_t);
    }
    std::cout << "user time = " 
              << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
              << "s\tsystem time = "
              << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
              << "s" << std::endl;
    
    // benchmark indexing with an indirect wrapped integer
    std::cout << "Testing indirect wrapped integer..." << std::endl;
    {
        Property<Integer> pint;
        times(&start_t);
        for (int i=0; i < NUM_TESTS; i++) {
            check<Integer, OtherInteger>(pint);
        }
        times(&stop_t);
    }
    std::cout << "user time = " 
              << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
              << "s\tsystem time = "
              << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
              << "s" << std::endl;
    
    return 0;
}

/*
  Local Variables:
  compile-command: "g++ -Wall -static -O3 \
  test_wrapped_int.C -o test_wrapped_int && ./test_wrapped_int"
  End:
*/
