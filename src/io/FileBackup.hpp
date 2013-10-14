#ifndef _IO_FILEBACKUP_HPP
#define _IO_FILEBACKUP_HPP
//#include <ctime>
#include <iostream>
//#include <sstream>
//#include <string>
//#include <boost/filesystem.hpp>

//namespace espressopp{
//  namespace io{
    class FileBackup {

      public:
        FileBackup(std::string file_name) { 
          std::cout << file_name << std::endl;
          /*

          //check if file with this name already exists
          if ( boost::filesystem::exists( file_name) )
          {
            //backup file
            //create backup filename %DATE%_FILE
            time_t t = time(0); // get time now
            struct tm * now = localtime (&t);
            std::stringstream new_file_name;
            new_file_name 
              << "2013"
              //<< std::to_string(now->tm_year + 1900)
              << "10"
              //<< std::to_string(now->tm_mon + 1)
              << (now->tm_mday)
              << (now->tm_hour)
              << (now->tm_min)
              << std::endl;
            std::cout << "Note: file " << file_name << " exists already. Moving " << file_name << " to " << new_file_name << std::endl;
            //make sure the new file does not exist already
          if ( boost::filesystem::exists( file_name) )
          {
            std::cout << "Warning: can not backup file " << file_name << " to " << new_file_name << ", because it exists already!" << std::endl;
            //die
          }
          else{
            //move FILE to %DATE%_FILE, where %DATE% is the current time and date
            boost::filesystem::copy_file(file_name, new_file_name.str());
            boost::filesystem::remove(file_name);
          }
        }
    */
    }
  };
//}
//}
#endif//_IO_FILEBACKUP_HPP
