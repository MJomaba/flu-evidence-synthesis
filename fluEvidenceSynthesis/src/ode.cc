#include "ode.h"

namespace ode {
    /*template<typename ODE_FUNC>
        Eigen::VectorXd step( Eigen::VectorXd &&y, ODE_FUNC &ode_func,
                double &step_size, const double start_time, 
                const double max_time )
        {
            step_size = std::min( step_size, max_time - start_time );
            return y + step_size*ode_func( y, start_time );
        }*/
};
