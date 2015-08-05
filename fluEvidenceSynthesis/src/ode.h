#ifndef FLU_ODE_H
#define FLU_ODE_H

//#include "rcppwrap.h"
#include<Eigen/Core>

namespace ode {
    template<typename ODE_FUNC>
        Eigen::VectorXd step( Eigen::VectorXd &&y, ODE_FUNC &ode_func,
                double &step_size, double &current_time, 
                const double max_time )
        {
            auto dt = std::min( step_size, max_time - current_time );
            y = y + dt*ode_func( y, current_time );
            current_time += dt;
            return y;
        }

    template<typename ODE_FUNC>
        Eigen::VectorXd rkf45_astep( Eigen::VectorXd &&y, 
                ODE_FUNC &ode_func,
                double &step_size, double &current_time, 
                const double max_time, const double tol = 1e-2 )
        {
            static const std::array<double, 4> tscale = 
                { 1.0/4, 3.0/8.0, 12.0/13.0, 1.0/2.0 };
            static const std::array<double, 2> k3scale = { 3.0/32, 9.0/32 };
            static const std::array<double, 3> k4scale = 
                { 1932.0/2197, 7200.0/2197,7296.0/2197 };
            static const std::array<double, 4> k5scale = 
                { 439.0/216, 8.0, 3680.0/513, 845.0/4104 };
            static const std::array<double, 5> k6scale = 
                { 8.0/27, 2.0, 3544.0/2565, 1859.0/4104,11.0/40 };
            static const std::array<double, 5> errscale = 
                { 1.0/360, 128.0/4275, 2197.0/75240, 1.0/50, 2.0/55 };
            static const std::array<double, 4> yscale =
                {25.0/216, 1408.0/2565, 2197.0/4101, 1.0/5};

            static Eigen::VectorXd k1, k2, k3, k4, k5, k6, err, y5;

            if (k1.size() != y.size())
            {
                k1 = Eigen::VectorXd( y.size() );
                k2 = k1;
                k3 = k1;
                k4 = k1;
                k5 = k1;
                k6 = k1;
                err = k1;
                y5 = k1;
            }

            auto max_step = max_time - current_time;

            auto dt = std::min( step_size, max_step );

            bool adapted = true;
            while( adapted )
            {

                k1 = dt*ode_func(y, current_time);
                k2 = dt*ode_func(
                        y + tscale[0]*k1, 
                        current_time + tscale[0]*dt);
                k3 = dt*ode_func(
                        y + k3scale[0]*k1+k3scale[1]*k2, 
                        current_time + tscale[1]*dt);
                k4 = dt*ode_func(
                        y + k4scale[0]*k1-k4scale[1]*k2+k4scale[2]*k3, 
                        current_time + tscale[2]*dt);
                k5 = dt*ode_func(
                        y + k5scale[0]*k1-k5scale[1]*k2+k5scale[2]*k3-
                        k5scale[3]*k4, 
                        current_time + dt);
                k6 = dt*ode_func(
                        y - k6scale[0]*k1+k6scale[1]*k2-k6scale[2]*k3+
                        k6scale[3]*k4-k6scale[4]*k5, 
                        current_time + tscale[3]*dt);

                err = ( errscale[0]*k1-errscale[1]*k3-errscale[2]*k4
                        +errscale[3]*k5+errscale[4]*k6 );

                auto s = std::min(
                        std::max(
                            0.84*pow(tol*dt/err.norm(),0.25), 
                            0.25)
                        , 5.0 );

                if ((s > 1.5 || s < 0.9) && dt < max_step)
                {
                    step_size = 0.9*s*dt; 
                    dt = std::min( step_size, max_step );
                } else {
                    adapted = false;
                }
            }

            // Solving some numerical problems
            if ( dt == max_step )
                current_time = max_time;
            else 
                current_time += dt;

            //Rcpp::cout << dt << ", " << step_size << ", " << current_time << std::endl;

            return y + yscale[0]*k1+yscale[1]*k3+yscale[2]*k4-
                yscale[3]*k5;
        }
};

#endif
