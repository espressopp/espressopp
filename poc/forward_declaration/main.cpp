#include "A.hpp"
#include "B.hpp"

int main() {
  B b;
  A a;

  a.b = &b;

  a.do_something();
}
