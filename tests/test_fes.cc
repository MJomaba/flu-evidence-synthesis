#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <iostream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>

#include "catch.hpp"

#include "model.hh"

namespace io = boost::iostreams;
namespace fs = boost::filesystem;

class FileWatcher {
    public:
        const std::string &_path;
        std::time_t lastWriteTime;

        FileWatcher( const std::string path )
            : _path( path )
        {
            if (!fs::exists( _path ))
                throw "No such file";
            std::cout << "Blaat" << std::endl;
            std::cout << _path << std::endl;
            lastWriteTime = fs::last_write_time( _path );
        }

        /// Has the file been modified since the watcher was created
        bool modified()
        {
            auto previousWriteTime = lastWriteTime;
            lastWriteTime = fs::last_write_time( _path );
            if (lastWriteTime > previousWriteTime)
                return true;
            return false;
        }
};

bool equal_file_contents( const std::string &path, 
        const std::string &other_path )
{
            io::mapped_file_source f1( path );
            io::mapped_file_source f2( other_path );

            if ( f1.size() > 0 
                    && f1.size() == f2.size() 
                    && std::equal( f1.data(), 
                        f1.data() + f1.size(), f2.data() ) 
               ) 
                return true;
            else
                return false;
}

TEST_CASE( "Run a short test run", "[full]" ) {
    auto fWatcher = FileWatcher( "../data/posterior.txt" );
    REQUIRE( system( "bin/flu-evidence-synthesis -d ../data/ --burn-in 1000 --chain-length 1000 --thinning 1" ) == 0 );
    REQUIRE( fWatcher.modified() );
    sleep(1);

    REQUIRE( equal_file_contents( 
                "../data/posterior.txt",
                "tests/short_posterior_sample.txt" ) );
}
