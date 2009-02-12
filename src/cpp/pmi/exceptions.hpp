#ifndef _PMI_EXCEPTIONS_HPP
#define _PMI_EXCEPTIONS_HPP
#include "pmi/types.hpp"

#include <exception>
#include <string>
#include <sstream>

using namespace std;

#define PMI_THROW_INTL_ERROR(output)		\
  {						\
    ostringstream ost; ost << output;		\
    LOG4ESPP_FATAL(logger, ost.str());		\
    throw InternalError(ost.str());		\
  }

#define PMI_THROW_USER_ERROR(output)		\
  {						\
    ostringstream ost; ost << output;		\
    LOG4ESPP_FATAL(logger, ost.str());		\
    throw UserError(ost.str());			\
  }

namespace pmi {
  class PMIException : public exception {};

  class UserError : public PMIException {
  public:
    string message;
    UserError(const string &_message);
    virtual const char* what() const throw();
    virtual ~UserError() throw ();
  };

  class InternalError : public PMIException {
  public:
    string message;
    InternalError(const string &_message);
    virtual const char* what() const throw();
    virtual ~InternalError() throw ();
  };

}
#endif
