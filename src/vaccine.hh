#ifndef FLU_VACCINE_HH
#define FLU_VACCINE_HH

#include <stdio.h>
#include <array>
#include <vector>

#include "io.hh"

namespace flu {
    namespace vaccine {

        /// Details of a vaccine programme
        struct vaccine_t {
            /// Efficacy for each year
            std::array<double, 7> efficacy_year;

            /// Calendar for the programme
            std::array<double, 2583> calendar;
        };

        /// Load a vaccine from the given file
        /// Assumes current line in FILE is where the file starts
        vaccine_t load_vaccine( FILE *file );

        /// Returns a vector with the whole vaccine programme
        /// 
        /// TODO: Do not rely on number of scenarios given in the beginning 
        std::vector<vaccine_t> load_vaccine_programme( 
                const std::string &path );
    };
};
#endif
