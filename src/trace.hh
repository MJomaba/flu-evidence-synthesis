/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Edwin van Leeuwen
  
  This file is part of chainmcmc.
  
  chainmcmc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  chainmcmc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with chainmcmc. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#ifndef TRACE_HH
#define TRACE_HH
#include<math.h>
#include<numeric>
#include<algorithm>
#include<vector>
#include<iostream>

#include <boost/lexical_cast.hpp>

namespace flu {
namespace trace {
	typedef std::vector<double> sample_t;

	struct sampleState {
		sample_t sample;
		double log_likelihood;
	};


	namespace details {
		double mean_v( const std::vector<double> &v );

		double var_v( const std::vector<double> &v );

		std::pair<double, double> confidence( const double interval,
				std::vector<double> v );

		double cov_v( const std::vector<double> &v,
				const std::vector<double> &w );	
	};

	/**
	 * \brief Reads a trace file
	 *
	 * Returns a vector with vectors of each parameter value
	 */
	std::vector<std::vector<double> > read_trace( 
			const std::string & fname );

	/**
	 * \brief Turns a string into a sample_t
	 */
	sample_t sample_from_string( const std::string & line );

	/**
	 * \brief Reads a trace file
	 *
	 * Returns a vector with samples
	 *
	 * \param tail Read the given number of samples from the end. Set to zero to read the whole file
	 */
	std::vector<std::vector<double> > read_trace_per_sample( 
			std::istream& infile, const size_t tail = 0 );

	/**
	 * \brief Read any new samples that have been appended to the trace file
	 */
	std::vector<std::vector<double> > follow_trace_per_sample( 
			std::istream& infile );

	/**
	 * \brief Returns means of the parameters in the trace
	 *
	 * Trace is assumed to be a vector of samples
	 */
	std::vector<double> means( const std::vector<sample_t> &samples );

	std::vector<double> variances_sample( 
			const std::vector<sample_t> & samples );
};
};

template<class T>
std::ostream& operator<<( std::ostream &out, const std::vector<T> &v ) {
	std::stringstream s; // Collect output in stringstream for thread safety
	bool first = true;
	/*std::cout << state.loglikelihood << ", " << temperature << 
		", " << state.current_parameter << ": ";*/
	for ( auto & el : v ) {
		if (first) {
			s << el;
			first = false;
		} else
			s	<< "\t" << el;
	}
	out << s.str();
	return out;
}
#endif
