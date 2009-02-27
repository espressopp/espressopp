#ifndef _ESUTIL_IDENTIFIER_HPP
#define _ESUTIL_IDENTIFIER_HPP

namespace esutil {
  /** class representing an identifier, essentially a nonmutable
      nonnegative integer. The identifier even allows 

      @tparam the class deriving from Identifier.
  */
  class Identifier {
  public:
    /// comparison
    bool operator==(const Identifier &o) const { return v == o.v; }
    /// comparison
    bool operator!=(const Identifier &o) const { return v != o.v; }

  protected:
    Identifier(size_t _v = 0): v(_v) {}

    size_t v;
  };
}

#endif
