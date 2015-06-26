#ifndef FLU_VACCINE_HH
#define FLU_VACCINE_HH

#include <stdio.h>
#include <array>
#include <vector>

namespace flu {
    namespace vaccine {
        /// Details of a vaccine programme
        struct vaccine_t {
            /// Efficacy for each year
            std::array<double, 7> efficacy_year;

            /// Calendar for the programme
            std::array<double, 2583> calendar;
        };
    };
};
#endif
