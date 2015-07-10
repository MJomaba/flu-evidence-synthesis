#ifndef FLU_VACCINE_HH
#define FLU_VACCINE_HH

#include <stdio.h>
#include <array>
#include <vector>
#include <boost/date_time.hpp>

#include<Eigen/Core>

namespace flu {
    namespace vaccine {
        /// Details of a vaccine programme
        struct vaccine_t {
            /// Efficacy for each year
            std::array<double, 7> efficacy_year;

            /// Calendar for the programme
            //std::array<double, 2583> calendar;
            Eigen::MatrixXd calendar;

            /**
             * \brief Holds the start date corresponding with each row in the calendar
             *
             * If length one longer than the number of rows, then the last
             * date is assumed to be the ending date of that vaccination rate
             * (set vaccination to zero after that date).
             *
             * If no dates are given at all, then the uk vaccination 
             * calendar is assumed (every day from october till February 
             * (not including Feb))
             */
            std::vector<boost::posix_time::ptime> dates;
        };
    };
};
#endif
