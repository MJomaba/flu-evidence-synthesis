#ifndef IO_HH 
#define IO_HH

#include <stdio.h>
#include <string>

namespace flu {
    FILE * read_file( const std::string path, const std::string filename );
    FILE * read_file( const std::string filename );
    FILE * write_file( const std::string filename );
    FILE * append_file( const std::string filename );

    void save_fgets( char * buffer, int size, FILE * file );
    void save_fscanf( FILE * stream, const char * format, ... );
};

#endif

