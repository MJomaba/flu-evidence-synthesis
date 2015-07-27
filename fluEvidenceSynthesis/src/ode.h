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

    template<typename ODE_FUNC>
        Eigen::VectorXd rkf45_astep( Eigen::VectorXd &&y, 
                ODE_FUNC &ode_func,
                double &step_size, const double start_time, 
                const double max_time, const double tol = 1e-2 )
        {
            Eigen::VectorXd k1, k2, k3, k4, k5, k6;


            step_size = std::min( step_size, max_time - start_time );

            bool adapted = true;
            while( adapted )
            {

                k1 = step_size*ode_func(y, start_time);
                k2 = step_size*ode_func(
                        y + 1.0/4.0*k1, 
                        start_time + 1.0/4.0*step_size);
                k3 = step_size*ode_func(
                        y + 3.0/32*k1+9.0/32*k2, 
                        start_time + 3.0/8.0*step_size);
                k4 = step_size*ode_func(
                        y + 1932.0/2197*k1-7200.0/2197*k2+7296.0/2197*k3, 
                        start_time + 12.0/13.0*step_size);
                k5 = step_size*ode_func(
                        y + 439.0/216*k1-8*k2+3680.0/513*k3-845.0/4104*k4, 
                        start_time + step_size);
                k6 = step_size*ode_func(
                        y - 8.0/27*k1+2.0*k2-3544.0/2565*k3+1859.0/4104*k4
                        -11.0/40*k5, 
                        start_time + 1.0/2.0*step_size);

                auto err = ( 1.0/360*k1-128.0/4275*k3-2197.0/75240*k4
                        +1.0/50*k5+2.0/55*k6 );

                auto s = std::min(
                        std::max(
                            0.84*pow(tol*step_size/err.norm(),0.25), 
                            0.25)
                        , 5.0 );

                if ((s > 2 || s < 0.5) && step_size < max_time-start_time)
                {
                    step_size = std::min( s*step_size, max_time - start_time );
                } else {
                    adapted = false;
                }
            }

            return y + 25.0/216*k1+1408.0/2565*k3+2197.0/4101*k4-
                1.0/5*k5;
        }
};

#endif
