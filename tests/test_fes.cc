#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <iostream>
#include <random>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "catch.hpp"

#include "model.hh"
#include "trace.hh"
#include "io.hh"
#include "vaccine.hh"
#include "proposal.hh"

namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace bu = boost::numeric::ublas;

using namespace flu;

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

bool approx_equal( double a, double b, double delta = 0 )
{
    if (std::abs( a-b ) <= delta)
        return true;
    return false;
}

bool equal_matrix_and_array( const bu::matrix<double> &a, double * b,
        double delta = 0 )
{
    bool equal = true;
    for( size_t i = 0; i < a.size1(); ++i )
    {
        for( size_t j = 0; j < a.size2(); ++j )
        {
            if ( !approx_equal(a(i,j), b[i*a.size2()+j], delta ) )
                equal = false;
        }
    }
    return equal;
}

bu::matrix<double> copy_matrix( double * a, size_t dim )
{
    bu::matrix<double> b(dim, dim);
    for( size_t i = 0; i < dim; ++i )
    {
        for( size_t j = 0; j < dim; ++j )
        {
            b(i,j) = a[i*dim+j];
        }
    }    
    return b;
}

TEST_CASE( "New method of loading vaccine data should give the same output",
        "[refactor]" )
{
    // Old way
    auto vacc_programme= read_file( "../data/vaccine_calendar.txt");
    char sbuffer[300];
    double vaccine_efficacy_year[7], vaccine_cal[2583], *VE_pro, *VCAL_pro;
    double *tab_cal[100], *tab_VE[100]; /*so far a maximum of 100 scenarios can be changed of course*/
    size_t i,j;
    size_t n_scenarii;

    save_fgets(sbuffer, 100, vacc_programme);
    save_fgets(sbuffer, 100, vacc_programme);
    sscanf(sbuffer,"%lu",&n_scenarii);

    save_fgets(sbuffer, 100, vacc_programme);
    save_fgets(sbuffer, 100, vacc_programme);
    save_fgets(sbuffer, 100, vacc_programme);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf",&vaccine_efficacy_year[0],&vaccine_efficacy_year[1],&vaccine_efficacy_year[2],&vaccine_efficacy_year[3],&vaccine_efficacy_year[4],&vaccine_efficacy_year[5],&vaccine_efficacy_year[6]);

    save_fgets(sbuffer, 50, vacc_programme);
    for(j=0;j<123;j++)
    {
        save_fgets(sbuffer, 300, vacc_programme);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &vaccine_cal[j*21],&vaccine_cal[j*21+1],&vaccine_cal[j*21+2],&vaccine_cal[j*21+3],&vaccine_cal[j*21+4],&vaccine_cal[j*21+5],&vaccine_cal[j*21+6],&vaccine_cal[j*21+7],&vaccine_cal[j*21+8],&vaccine_cal[j*21+9],&vaccine_cal[j*21+10],&vaccine_cal[j*21+11],&vaccine_cal[j*21+12],&vaccine_cal[j*21+13],&vaccine_cal[j*21+14],&vaccine_cal[j*21+15],&vaccine_cal[j*21+16],&vaccine_cal[j*21+17],&vaccine_cal[j*21+18],&vaccine_cal[j*21+19],&vaccine_cal[j*21+20]);
    }

    tab_cal[0]=vaccine_cal;
    tab_VE[0]=vaccine_efficacy_year;

    for(i=0;i<n_scenarii;i++)
    {
        VE_pro=(double *) malloc(7 * sizeof (double));
        save_fgets(sbuffer, 100, vacc_programme);
        save_fgets(sbuffer, 100, vacc_programme);
        save_fgets(sbuffer, 100, vacc_programme);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf",&VE_pro[0],&VE_pro[1],&VE_pro[2],&VE_pro[3],&VE_pro[4],&VE_pro[5],&VE_pro[6]);
        tab_VE[1+i]=VE_pro;

        VCAL_pro=(double *) malloc(2583 * sizeof (double));
        save_fgets(sbuffer, 100, vacc_programme);
        for(j=0;j<123;j++)
        {
            save_fgets(sbuffer, 300, vacc_programme);
            sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &VCAL_pro[j*21],&VCAL_pro[j*21+1],&VCAL_pro[j*21+2],&VCAL_pro[j*21+3],&VCAL_pro[j*21+4],&VCAL_pro[j*21+5],&VCAL_pro[j*21+6],&VCAL_pro[j*21+7],&VCAL_pro[j*21+8],&VCAL_pro[j*21+9],&VCAL_pro[j*21+10],&VCAL_pro[j*21+11],&VCAL_pro[j*21+12],&VCAL_pro[j*21+13],&VCAL_pro[j*21+14],&VCAL_pro[j*21+15],&VCAL_pro[j*21+16],&VCAL_pro[j*21+17],&VCAL_pro[j*21+18],&VCAL_pro[j*21+19],&VCAL_pro[j*21+20]);
        }
        tab_cal[1+i]=VCAL_pro;
    }

    auto vaccine_programme = vaccine::load_vaccine_programme( "../data/vaccine_calendar.txt" );

    for( size_t i = 0; i < vaccine_programme.size(); ++i )
    {
        for ( size_t j = 0; j < vaccine_programme[i].efficacy_year.size();
                ++j )
        {
            INFO( "Efficacy i,j: " << i << "," << j );
            REQUIRE( vaccine_programme[i].efficacy_year[j] ==
                tab_VE[i][j] );
        }
        for ( size_t j = 0; j < vaccine_programme[i].calendar.size();
                ++j )
        {
            INFO( "Calendar i,j: " << i << "," << j );
            REQUIRE( vaccine_programme[i].calendar[j] ==
                tab_cal[i][j] );
        }    
    }
}

TEST_CASE( "New cholesky should give the same output", "[refactor]" )
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution( 0, 10 );
    double mat[25];
    double mat2[25];
    for (size_t i=0; i<5; ++i) {
        for (size_t j=0; j<5; ++j) {
            if (j<=i)
                mat[i*5+j] = distribution(generator);
            else
                mat[i*5+j] = mat[j*5+i];
        }
    }
    auto bmat = copy_matrix( mat, 5 );
    REQUIRE( equal_matrix_and_array( bmat, mat ) );

    flu::cholevsky( mat, mat2, 5 );
    auto bmat2 = flu::proposal::cholesky_factorization( bmat );
    REQUIRE( equal_matrix_and_array( bmat2, mat2, 1e-100 ) );
}

TEST_CASE( "Run a short test run", "[hide]" ) 
{
    // Create array with all names to be watched
    std::vector<std::string> files;
    files.push_back( "posterior.txt" );
    files.push_back( "final.log" );

    // Add sample files to check
    for( size_t i = 1010; i < 2001; i+= 10 )
    {
        files.push_back( "samples/z_hyper" 
                + boost::lexical_cast<std::string>( i ) + ".stm" );
    }

    // Start with clean slate
    for( auto &fname : files )
    {
        fs::remove( "../data/" + fname );
    }


    REQUIRE( system( "bin/flu-evidence-synthesis -d ../data/ --burn-in 1000 --chain-length 1000 --thinning 1" ) == 0 );

    for( auto &fname : files )
    {
        REQUIRE( fs::exists("../data/" + fname) );

        auto full_path = "../data/" + fname;
        auto test_path = "./tests/test_data/" + fname;

        INFO( "Files not equal: " << full_path << " " << test_path );

        REQUIRE( equal_file_contents( 
                    full_path, test_path ) );
    }
}

TEST_CASE( "Run inference", "[full,inference]" )
{
    // Create array with all names to be watched
    std::vector<std::string> files;

    // Setup the scenarii for watching
    std::string pre = "scenarii/Scenario_";
    std::string post = "_final_size.txt";
    files.push_back( pre + "vaccination" + post );
    files.push_back( pre + "no_vaccination" + post );
    for ( size_t i = 0; i<35; ++i )
        files.push_back( pre + boost::lexical_cast<std::string>( i )
                + post );

    // Start with clean slate
    for( auto &fname : files )
    {
        fs::remove( "../data/" + fname );
    }


    REQUIRE( system( "bin/inference -d ../data/" ) == 0 );

    for( auto &fname : files )
    {
        REQUIRE( fs::exists("../data/" + fname) );

        auto full_path = "../data/" + fname;
        auto test_path = "./tests/test_data/" + fname;

        /*INFO( "Files not equal: " << full_path << " " << test_path );
        std::ifstream infile;
        infile.open( full_path );
        auto new_scenario = flu::trace::read_trace_per_sample(infile);
        std::ifstream infile2;
        infile2.open( test_path );
        auto orig_scenario = flu::trace::read_trace_per_sample(infile2);

        for( size_t i = 0; i <new_scenario.size(); ++i )
        {
            for ( size_t j = 0; j < new_scenario[i].size(); ++j )
            {
                if (new_scenario[i][j] == 0)
                    REQUIRE( orig_scenario[i][j] == 0 );
                else {
                    double rel_diff = std::abs(
                            (orig_scenario[i][j]-new_scenario[i][j])
                            /new_scenario[i][j]);
                    std::cout << rel_diff << std::endl;
                    REQUIRE( rel_diff <= 0.0001 );
                }
            }
        }*/

        INFO( "Files not equal: " << full_path << " " << test_path );
        REQUIRE( equal_file_contents( 
                    full_path, test_path ) );
    }
}

TEST_CASE( "Posterior file should be readable as a trace" "[trace]" )
{
    std::string fName = "./tests/test_data/posterior.txt";

    std::ifstream infile;
    infile.open( fName );

    std::vector<std::vector<double> > samples 
        = flu::trace::read_trace_per_sample( infile, 1 );

    REQUIRE( samples.size() == 1000 );
    REQUIRE( samples[0].size() == 12 );

    infile.close();
}

