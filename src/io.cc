#include "io.hh"

#include <cstdarg>
#include<fstream>
#include<iostream>

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

    FILE * read_file( const std::string filename )
    {
        boost::filesystem::path path = filename;
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

	void move_last_lines( std::istream& fin, const size_t n ) {
		fin.seekg( 0,std::ios_base::end);
		if (n!=0) { 
			size_t no_found = 0;
			while (no_found < n) {
				fin.unget(); fin.unget(); // Move backwards two places

				char ch;
				fin.get(ch); // Get current char and move forward a place
				if ( (int)fin.tellg() <=1) { // We're at beginning of file
					fin.seekg(0);
					no_found = n; // Interrupt the loop
				} else {
					if (ch == '\n') {
						++no_found;
					}	
				}
			}
		}
	}

	std::vector<std::string> get_lines_and_move( std::istream& infile ) {
		std::vector<std::string> lines;
		std::string line;
		while (std::getline(infile, line)) {
			lines.push_back( line );
		}
		infile.clear(); // unset eofbit (and/or failbit)
		return lines;
	}

	std::vector<std::string> tail( const std::string & fname, 
			const size_t n ) {
		std::vector<std::string> lines;
		std::ifstream fin( fname );

		// Special rule that if n==0 we return whole file
		if (n!=0) {
			move_last_lines( fin, n );
		}
		
		return get_lines_and_move( fin );
	}
};
