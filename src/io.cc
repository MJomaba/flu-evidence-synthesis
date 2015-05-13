#include "io.hh"

#include <cstdarg>

#include "boost/filesystem.hpp"

namespace flu {

    FILE * read_file( const std::string pathname, const std::string filename )
    {
        boost::filesystem::path path = pathname;
        path /= filename;
        if (!boost::filesystem::exists( path )) {
            std::cerr << "File does not exist: " << path << std::endl;
        }
        return fopen(path.c_str(),"r");
    }

    FILE * write_file( const std::string filename )
    {
        boost::filesystem::path filepath = filename;
        boost::filesystem::path path = filename;
        // Create directory if it doesn't exist
        // TODO: There is probably a function to get the path without using remove_filename()
        // If so then the copy is not needed anymore
        boost::filesystem::create_directory( path.remove_filename() );
        FILE * file = fopen( filepath.c_str(), "w+t" );
        return file;
    }

    FILE * append_file( const std::string filename )
    {
        boost::filesystem::path filepath = filename;
        boost::filesystem::path path = filename;
        // Create directory if it doesn't exist
        // TODO: There is probably a function to get the path without using remove_filename()
        // If so then the copy is not needed anymore
        boost::filesystem::create_directory( path.remove_filename() );
        if (!boost::filesystem::exists( filepath )) {
            return fopen(filepath.c_str(), "w+t");
        }
        return fopen(filepath.c_str(), "a");
    }

    //! Check for return value of fgets 
    void save_fgets( char * buffer, int size, FILE * file ) 
    {
        auto res = fgets( buffer, size, file );
        if (!res)
            throw "Could not read from file";
    }

    //! Check for return value of fscanf 
    void save_fscanf( FILE * stream, const char * format, ... ) 
    {
        va_list args;
        va_start( args, format );
        auto res = vfscanf( stream, format, args );
        va_end( args );
        if (!res)
            throw "Could not read from file";
    }
};
