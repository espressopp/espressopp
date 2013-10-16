#ifndef _IO_FILEBACKUP_HPP
#define _IO_FILEBACKUP_HPP
#include <iostream>
#include <sstream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include "System.hpp"


    class FileBackup {

      public:
        FileBackup(std::string file_name) { 

          //check if file with this name already exists
          if ( boost::filesystem::exists( file_name) )
          {
            //backup file
            //create backup filename %DATE%_FILE
            boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
            std::stringstream ssnew_file_name;
            ssnew_file_name 
              << now.date().year()
              << now.date().month()
              << now.date().day()
              << now.time_of_day().hours()
              << now.time_of_day().minutes()
              << "_" << file_name; // this is a bit ugly since there are no leading zeros
            std::string new_file_name(ssnew_file_name.str()); 
            std::cout << "Note: file " << file_name << " exists already. Moving " << file_name << " to " << new_file_name << std::endl;
            //make sure the new file does not exist already
          if ( boost::filesystem::exists( new_file_name) )
          {
            //die
            std::stringstream msg;
            msg << "Warning: can not backup file " << file_name << " to " << new_file_name << ", because it exists already!";
            throw std::runtime_error(msg.str());
          }
          else{
            //move FILE to %DATE%_FILE, where %DATE% is the current time and date
            boost::filesystem::rename(file_name, new_file_name);
          }
        }
    }
  };

#endif//_IO_FILEBACKUP_HPP
