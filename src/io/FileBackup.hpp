/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

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
