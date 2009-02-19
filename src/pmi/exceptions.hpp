#ifndef _PMI_EXCEPTIONS_HPP
#define _PMI_EXCEPTIONS_HPP
#include "pmi/types.hpp"

#include <exception>
#include <string>
#include <sstream>

#define PMI_THROW_INTL_ERROR(output)		\
  {						\
    std::ostringstream ost; ost << output;	\
    LOG4ESPP_FATAL(pmi::logger, ost.str());	\
    throw InternalError(ost.str());		\
  }

#define PMI_THROW_USER_ERROR(output)		\
  {						\
    std::ostringstream ost; ost << output;	\
    LOG4ESPP_FATAL(pmi::logger, ost.str());	\
    throw UserError(ost.str());			\
  }

namespace pmi {
  class PMIException : public std::exception {};

  class UserError : public PMIException {
  public:
    std::string message;
    UserError(const std::string &_message);
    virtual const char* what() const throw();
    virtual ~UserError() throw ();
  };

  class InternalError : public PMIException {
  public:
    std::string message;
    InternalError(const std::string &_message);
    virtual const char* what() const throw();
    virtual ~InternalError() throw ();
  };

}
#endif
