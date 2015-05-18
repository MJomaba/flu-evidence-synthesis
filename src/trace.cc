/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Edwin van Leeuwen
  
  This file is part of chainmcmc.
  
  chainmcmc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Chainmcmc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Gillespie. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include<fstream>

#include <boost/algorithm/string.hpp>

#include "trace.hh"
#include "io.hh"
namespace flu {
	namespace trace {
		namespace details {
			double mean_v( const std::vector<double> &v ) {
				double sum = std::accumulate( v.begin(), v.end(), 0.0, 
						[]( const double a, const double b ) { return a+b; } );
				return sum/v.size();
			}

			double var_v( const std::vector<double> &v ) {
				double mean = mean_v( v );
				double var_sum = std::accumulate( v.begin(), v.end(), 0.0, 
						[&mean]( const double a, const double b ) { 
						return a+pow(b-mean,2); } );
				return var_sum/v.size();
			}

			std::pair<double, double> confidence( const double interval,
					std::vector<double> v ) {
				std::pair<double, double> result;
				std::sort( v.begin(), v.end() );
				double alpha = (1.0-interval)/2.0;
				result.first = v[ceil(alpha*v.size())];
				result.second = v[floor((1.0-alpha)*v.size())];
				return result;
			}

			double cov_v( const std::vector<double> &v,
					const std::vector<double> &w ) {
				if ( v.size() != w.size() ) {
					std::cerr << "Vectors need to be the same size" << std::endl;
					throw;
				}
				double meanv = mean_v( v );
				double meanw = mean_v( w );
				double cov = 0;
				for ( size_t i = 0; i < v.size(); ++i ) {
					cov += (v[i] - meanv)*(w[i]-meanw);
				}
				return cov/v.size();
			}
		};

		std::vector<double> sample_from_string( const std::string & line ) {
			std::vector<std::string> strs;
			boost::split( strs, line, boost::is_any_of( " " ) );

			std::vector<double> pars;
			for ( auto & str : strs )
				pars.push_back( boost::lexical_cast<double>( str ) );
			return pars;
		}

		std::vector<std::vector<double> > read_trace_per_sample( 
				std::istream& infile, size_t skip, const size_t tail ) {
			if (tail != 0) {
				move_last_lines( infile, tail );
			}
			return follow_trace_per_sample( infile, skip );
		}

		std::vector<std::vector<double> > read_trace( 
				const std::string & fname ) {
			std::vector<std::vector<double> > trace;
			std::ifstream infile;
		infile.open( fname );
		auto totranspose = read_trace_per_sample( infile, 0 );
		// Now transpose
		for ( auto & sample : totranspose ) {
			for ( size_t i = 0; i < sample.size(); ++i ) {
				if (trace.size()<i)
					trace.push_back( std::vector<double>() );
				trace[i].push_back( sample[i] );
			}
		}
		return trace;
	}

	std::vector<std::vector<double> > follow_trace_per_sample(
		std::istream& infile, size_t skip ) {
		std::vector<std::vector<double> > samples;
        auto lines = get_lines_and_move( infile );
		for ( auto & line : lines ) {
            if (skip <= 0)
                samples.push_back( sample_from_string( line ) );
            else
                --skip;
		}
		return samples;
	}

	std::vector<double> means( const std::vector<sample_t> & samples ) {
		std::vector<double> means;
		if (samples.size()>0) {
			means = std::vector<double>( samples[0].size() );
			for ( auto & mean : means )
				mean = 0;
			for ( auto & sample : samples ) {
				for ( size_t i = 0; i < sample.size(); ++i ) {
					means[i] += sample[i];
				}
			}
		}
		for ( auto & mean : means )
			mean /= samples.size();
		return means;
	}

	std::vector<double> variances_sample( 
			const std::vector<sample_t> & samples ) {
		std::vector<double> ms = means( samples );
		std::vector<double> vs( ms.size() );

		if (samples.size()>0) {
			for ( auto & v : vs )
				v = 0;
			for ( auto & sample : samples ) {
				for ( size_t i = 0; i < sample.size(); ++i ) {
					vs[i] += pow(ms[i] - sample[i],2);
				}
			}
		}
		for ( auto & v : vs )
			v /= samples.size();
		return vs;
	}	
	};
};
