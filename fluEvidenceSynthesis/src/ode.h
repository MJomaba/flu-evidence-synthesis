#ifndef FLU_ODE_H
#define FLU_ODE_H

#include<Eigen/Core>

namespace ode {
    template<typename ODE_FUNC>
        Eigen::VectorXd step( Eigen::VectorXd &&y, ODE_FUNC &ode_func,
                double &step_size, const double start_time, 
                const double max_time )
        {
            step_size = std::min( step_size, max_time - start_time );
            return y + step_size*ode_func( y, start_time );
        }

    /*template<typename Func>
        Eigen::VectorXd solve( Eigen::VectorXd &&y, Func &ode_func,
                const boost::posix_time::ptime &start_time,
                const boost::posix_time::ptime &end_time );*/
};

#endif
