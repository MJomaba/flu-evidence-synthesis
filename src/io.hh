#ifndef IO_HH 
#define IO_HH

#include <stdio.h>
#include <string>
#include<vector>

namespace flu {
    FILE * read_file( const std::string path, const std::string filename );
    FILE * read_file( const std::string filename );
    FILE * write_file( const std::string filename );
    FILE * append_file( const std::string filename );

    void save_fgets( char * buffer, int size, FILE * file );
    void save_fscanf( FILE * stream, const char * format, ... );

	/**
	 * \brief Move to the beginning of the last n lines
	 *
	 * Implementation detail: cannot return an ifstream, not sure if move is still not
	 * implemented in c++11 or if I am missing how exactly to do it.
	 */
	void move_last_lines( std::istream& infile, const size_t n );

	/**
	 * \brief Get lines starting from current position and move position to the end
	 */
	std::vector<std::string> get_lines_and_move( std::istream& infile );

	/**
	 * \brief Read last lines from the given file
	 */
	std::vector<std::string> tail( const std::string & fname, 
			const size_t no );
};

#endif

