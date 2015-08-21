#include "model.h"

#include "ode.h"

namespace flu 
{
    boost::posix_time::ptime getTimeFromWeekYear( int week, int year )
    {
        namespace bt = boost::posix_time;
        namespace bg = boost::gregorian;
        // Week 1 is the first week that ends in this year
        auto firstThursday = bg::first_day_of_the_week_in_month( 
                bg::Thursday, bg::Jan );
        auto dateThursday = firstThursday.get_date( year );
        auto current_time = bt::ptime(dateThursday) - bt::hours(24*3);
        current_time += bt::hours(24*7*(week-1));
        return current_time;
    }

    enum seir_type_t { S = 0, E1 = 1, E2 = 2, I1 = 3, I2 = 4, R = 5 };
    const std::vector<seir_type_t> seir_types = 
        { S, E1, E2, I1, I2, R };

    enum group_type_t { LOW = 0, HIGH = 1, PREG = 2,
        VACC_LOW = 3, VACC_HIGH = 4, VACC_PREG = 5 };
    const std::vector<group_type_t> group_types = { LOW, HIGH, PREG,
        VACC_LOW, VACC_HIGH, VACC_PREG };

    inline size_t ode_id( const size_t nag, const group_type_t gt, 
            const seir_type_t st )
    {
        return gt*nag*seir_types.size()
        + st*nag;
    }

    inline size_t ode_id( const size_t nag, const group_type_t gt, 
            const seir_type_t st, 
            const size_t i )
    {
        return gt*nag*seir_types.size()
        + st*nag + i;
    }

    inline Eigen::VectorXd flu_ode( Eigen::VectorXd &deltas,
            const Eigen::VectorXd &densities,
            const std::vector<double> &Npop,
            const Eigen::MatrixXd &vaccine_rates, // If empty, rate of zero is assumed
            const std::array<double,7> &vaccine_efficacy_year,
            const Eigen::MatrixXd &transmission_regular,
            double a1, double a2, double g1, double g2 )
    {
        const size_t nag = transmission_regular.cols();
        for(size_t i=0;i<nag;i++)
        {
            /*rate of depletion of susceptible*/
            deltas[ode_id(nag,VACC_LOW,S,i)]=0;
            for(size_t j=0;j<nag;j++)
                deltas[ode_id(nag,VACC_LOW,S,i)]+=transmission_regular(i,j)*(densities[ode_id(nag,VACC_LOW,I1,j)]+densities[ode_id(nag,VACC_LOW,I2,j)]+densities[ode_id(nag,VACC_HIGH,I1,j)]+densities[ode_id(nag,VACC_HIGH,I2,j)]+densities[ode_id(nag,VACC_PREG,I1,j)]+densities[ode_id(nag,VACC_PREG,I2,j)]+densities[ode_id(nag,LOW,I1,j)]+densities[ode_id(nag,LOW,I2,j)]+densities[ode_id(nag,HIGH,I1,j)]+densities[ode_id(nag,HIGH,I2,j)]+densities[ode_id(nag,PREG,I1,j)]+densities[ode_id(nag,PREG,I2,j)]);

            deltas[ode_id(nag,VACC_HIGH,S,i)]=deltas[ode_id(nag,VACC_LOW,S,i)];
            deltas[ode_id(nag,VACC_PREG,S,i)]=deltas[ode_id(nag,VACC_LOW,S,i)];
            deltas[ode_id(nag,LOW,S,i)]=deltas[ode_id(nag,VACC_LOW,S,i)];
            deltas[ode_id(nag,HIGH,S,i)]=deltas[ode_id(nag,VACC_LOW,S,i)];
            deltas[ode_id(nag,PREG,S,i)]=deltas[ode_id(nag,VACC_LOW,S,i)];

            deltas[ode_id(nag,VACC_LOW,S,i)]*=-densities[ode_id(nag,VACC_LOW,S,i)];
            deltas[ode_id(nag,VACC_HIGH,S,i)]*=-densities[ode_id(nag,VACC_HIGH,S,i)];
            deltas[ode_id(nag,VACC_PREG,S,i)]*=-densities[ode_id(nag,VACC_PREG,S,i)];
            deltas[ode_id(nag,LOW,S,i)]*=-densities[ode_id(nag,LOW,S,i)];
            deltas[ode_id(nag,HIGH,S,i)]*=-densities[ode_id(nag,HIGH,S,i)];
            deltas[ode_id(nag,PREG,S,i)]*=-densities[ode_id(nag,PREG,S,i)];
        }

        /*rate of passing between states of infection*/
        for ( auto &gt : group_types)
        {
            deltas.segment(ode_id(nag,gt,E1),nag)=-deltas.segment(ode_id(nag,gt,S),nag)-a1*densities.segment(ode_id(nag,gt,E1),nag);
            deltas.segment(ode_id(nag,gt,E2),nag)=a1*densities.segment(ode_id(nag,gt,E1),nag)-a2*densities.segment(ode_id(nag,gt,E2),nag);

            deltas.segment(ode_id(nag,gt,I1),nag)=a2*densities.segment(ode_id(nag,gt,E2),nag)-g1*densities.segment(ode_id(nag,gt,I1),nag);
            deltas.segment(ode_id(nag,gt,I2),nag)=g1*densities.segment(ode_id(nag,gt,I1),nag)-g2*densities.segment(ode_id(nag,gt,I2),nag);
            deltas.segment(ode_id(nag,gt,R),nag)=g2*densities.segment(ode_id(nag,gt,I2),nag);
        }

        /*Vaccine bit*/

        if ( vaccine_rates.size() > 0 )
        {
            for(size_t i=0;i<nag;i++)
            {
                double vacc_prov=Npop[i]*vaccine_rates(i)/(densities[ode_id(nag,LOW,S,i)]+densities[ode_id(nag,LOW,E1,i)]+densities[ode_id(nag,LOW,E2,i)]+densities[ode_id(nag,LOW,I1,i)]+densities[ode_id(nag,LOW,I2,i)]+densities[ode_id(nag,LOW,R,i)]);
                /*surv[i]+=vaccination_calendar[cal_time*21+i];*/
                double vacc_prov_r=Npop[i+nag]*vaccine_rates(i+nag)/(densities[ode_id(nag,HIGH,S,i)]+densities[ode_id(nag,HIGH,E1,i)]+densities[ode_id(nag,HIGH,E2,i)]+densities[ode_id(nag,HIGH,I1,i)]+densities[ode_id(nag,HIGH,I2,i)]+densities[ode_id(nag,HIGH,R,i)]);
                double vacc_prov_p=0; /*Npop[i+2*nag]*vaccination_calendar[cal_time*21+i+2*nag]/(densities[ode_id(nag,PREG,S,i)]+densities[ode_id(nag,PREG,E1,i)]+densities[ode_id(nag,PREG,E2,i)]+densities[ode_id(nag,PREG,I1,i)]+densities[ode_id(nag,PREG,I2,i)]+densities[ode_id(nag,PREG,R,i)]);*/

                deltas[ode_id(nag,VACC_LOW,S,i)]+=densities[ode_id(nag,LOW,S,i)]*vacc_prov*(1-vaccine_efficacy_year[i]);
                deltas[ode_id(nag,VACC_HIGH,S,i)]+=densities[ode_id(nag,HIGH,S,i)]*vacc_prov_r*(1-vaccine_efficacy_year[i]);
                deltas[ode_id(nag,VACC_PREG,S,i)]+=densities[ode_id(nag,PREG,S,i)]*vacc_prov_p*(1-vaccine_efficacy_year[i]);
                deltas[ode_id(nag,LOW,S,i)]-=densities[ode_id(nag,LOW,S,i)]*vacc_prov;
                deltas[ode_id(nag,HIGH,S,i)]-=densities[ode_id(nag,HIGH,S,i)]*vacc_prov_r;
                deltas[ode_id(nag,PREG,S,i)]-=densities[ode_id(nag,PREG,S,i)]*vacc_prov_p;

                deltas[ode_id(nag,VACC_LOW,E1,i)]+=densities[ode_id(nag,LOW,E1,i)]*vacc_prov;
                deltas[ode_id(nag,VACC_HIGH,E1,i)]+=densities[ode_id(nag,HIGH,E1,i)]*vacc_prov_r;
                deltas[ode_id(nag,VACC_PREG,E1,i)]+=densities[ode_id(nag,PREG,E1,i)]*vacc_prov_p;
                deltas[ode_id(nag,LOW,E1,i)]-=densities[ode_id(nag,LOW,E1,i)]*vacc_prov;
                deltas[ode_id(nag,HIGH,E1,i)]-=densities[ode_id(nag,HIGH,E1,i)]*vacc_prov_r;
                deltas[ode_id(nag,PREG,E1,i)]-=densities[ode_id(nag,PREG,E1,i)]*vacc_prov_p;

                deltas[ode_id(nag,VACC_LOW,E2,i)]+=densities[ode_id(nag,LOW,E2,i)]*vacc_prov;
                deltas[ode_id(nag,VACC_HIGH,E2,i)]+=densities[ode_id(nag,HIGH,E2,i)]*vacc_prov_r;
                deltas[ode_id(nag,VACC_PREG,E2,i)]+=densities[ode_id(nag,PREG,E2,i)]*vacc_prov_p;
                deltas[ode_id(nag,LOW,E2,i)]-=densities[ode_id(nag,LOW,E2,i)]*vacc_prov;
                deltas[ode_id(nag,HIGH,E2,i)]-=densities[ode_id(nag,HIGH,E2,i)]*vacc_prov_r;
                deltas[ode_id(nag,PREG,E2,i)]-=densities[ode_id(nag,PREG,E2,i)]*vacc_prov_p;

                deltas[ode_id(nag,VACC_LOW,I1,i)]+=densities[ode_id(nag,LOW,I1,i)]*vacc_prov;
                deltas[ode_id(nag,VACC_HIGH,I1,i)]+=densities[ode_id(nag,HIGH,I1,i)]*vacc_prov_r;
                deltas[ode_id(nag,VACC_PREG,I1,i)]+=densities[ode_id(nag,PREG,I1,i)]*vacc_prov_p;
                deltas[ode_id(nag,LOW,I1,i)]-=densities[ode_id(nag,LOW,I1,i)]*vacc_prov;
                deltas[ode_id(nag,HIGH,I1,i)]-=densities[ode_id(nag,HIGH,I1,i)]*vacc_prov_r;
                deltas[ode_id(nag,PREG,I1,i)]-=densities[ode_id(nag,PREG,I1,i)]*vacc_prov_p;

                deltas[ode_id(nag,VACC_LOW,I2,i)]+=densities[ode_id(nag,LOW,I2,i)]*vacc_prov;
                deltas[ode_id(nag,VACC_HIGH,I2,i)]+=densities[ode_id(nag,HIGH,I2,i)]*vacc_prov_r;
                deltas[ode_id(nag,VACC_PREG,I2,i)]+=densities[ode_id(nag,PREG,I2,i)]*vacc_prov_p;
                deltas[ode_id(nag,LOW,I2,i)]-=densities[ode_id(nag,LOW,I2,i)]*vacc_prov;
                deltas[ode_id(nag,HIGH,I2,i)]-=densities[ode_id(nag,HIGH,I2,i)]*vacc_prov_r;
                deltas[ode_id(nag,PREG,I2,i)]-=densities[ode_id(nag,PREG,I2,i)]*vacc_prov_p;

                deltas[ode_id(nag,VACC_LOW,R,i)]+=densities[ode_id(nag,LOW,R,i)]*vacc_prov+densities[ode_id(nag,LOW,S,i)]*vacc_prov*vaccine_efficacy_year[i];
                deltas[ode_id(nag,VACC_HIGH,R,i)]+=densities[ode_id(nag,HIGH,R,i)]*vacc_prov_r+densities[ode_id(nag,HIGH,S,i)]*vacc_prov_r*vaccine_efficacy_year[i];
                deltas[ode_id(nag,VACC_PREG,R,i)]+=densities[ode_id(nag,PREG,R,i)]*vacc_prov_p+densities[ode_id(nag,PREG,S,i)]*vacc_prov_p*vaccine_efficacy_year[i];
                deltas[ode_id(nag,LOW,R,i)]-=densities[ode_id(nag,LOW,R,i)]*vacc_prov;
                deltas[ode_id(nag,HIGH,R,i)]-=densities[ode_id(nag,HIGH,R,i)]*vacc_prov_r;
                deltas[ode_id(nag,PREG,R,i)]-=densities[ode_id(nag,PREG,R,i)]*vacc_prov_p;
            }
        }
        return deltas;
    }

    inline Eigen::VectorXd new_cases( 
            Eigen::VectorXd &densities,
            const boost::posix_time::ptime &start_time,
            const boost::posix_time::ptime &end_time, 
            boost::posix_time::time_duration &dt,
            const std::vector<double> &Npop,
            const Eigen::MatrixXd &vaccine_rates, // If empty, rate of zero is assumed
            const std::array<double,7> &vaccine_efficacy_year,
            const Eigen::MatrixXd &transmission_regular,
            double a1, double a2, double g1, double g2
            )
    {
        namespace bt = boost::posix_time;

        double h_step = dt.hours()/24.0;

        const size_t nag = transmission_regular.cols();
        Eigen::VectorXd results = Eigen::VectorXd::Zero(nag*3);

        static Eigen::VectorXd deltas( nag*group_types.size()*
                seir_types.size() );

        auto time_left = (end_time-start_time).hours()/24.0;

        auto ode_func = [&]( const Eigen::VectorXd &y, const double dummy )
        {
            return flu_ode( deltas, y, 
                    Npop, vaccine_rates, vaccine_efficacy_year,
                    transmission_regular, a1, a2, g1, g2 );
        };

        auto t = 0.0;
        auto prev_t = t;

        while (t < time_left)
        {
            prev_t = t;
            /*densities = ode::rkf45_astep( std::move(densities), ode_func,
                        h_step, t, time_left, 5 );*/
            densities = ode::step( std::move(densities), ode_func,
                        h_step, t, time_left );
            //Rcpp::Rcout << h_step << std::endl;

            results.block( 0, 0, nag, 1 ) += a2*(densities.segment(ode_id(nag,VACC_LOW,E2),nag)+densities.segment(ode_id(nag,LOW,E2),nag))*(t-prev_t);
            results.block( nag, 0, nag, 1 ) += a2*(densities.segment(ode_id(nag,VACC_HIGH,E2),nag)+densities.segment(ode_id(nag,HIGH,E2),nag))*(t-prev_t);
            results.block( 2*nag, 0, nag, 1 ) += a2*(densities.segment(ode_id(nag,VACC_PREG,E2),nag)+densities.segment(ode_id(nag,PREG,E2),nag))*(t-prev_t);
        }
        return results;
    }

    cases_t one_year_SEIR_with_vaccination(
            const std::vector<double> &Npop, double * seeding_infectious, 
            const double tlatent, const double tinfectious, 
            const std::vector<double> &s_profile, 
            const Eigen::MatrixXd &contact_regular, double q,
            const vaccine::vaccine_t &vaccine_programme,
            size_t minimal_resolution )
    {
        namespace bt = boost::posix_time;

        const size_t nag = contact_regular.rows(); // No. of age groups

        Eigen::VectorXd densities = Eigen::VectorXd::Zero( 
                nag*group_types.size()*
                seir_types.size() );

        double a1, a2, g1, g2 /*, surv[7]={0,0,0,0,0,0,0}*/;

        // We start at week 35. Week 1 is the first week that ends in this year
        auto current_time = getTimeFromWeekYear( 35, 1970 );
        if (vaccine_programme.dates.size()!=0)
            current_time = getTimeFromWeekYear( 35, 
                    vaccine_programme.dates[0].date().year() );

        auto start_time = current_time;
        auto end_time = current_time + bt::hours(364*24);
        size_t step_count = 0;

        a1=2/tlatent;
        a2=a1;
        g1=2/tinfectious;
        g2=g1;

        int date_id = -1;

        /*initialisation, transmission matrix*/
        Eigen::MatrixXd transmission_regular(contact_regular);
        for(int i=0;i<transmission_regular.rows();i++)
        {
            for(int j=0;j<transmission_regular.cols();j++) {
                transmission_regular(i,j)*=q*s_profile[i];
            }
        }

        /*initialisation, densities.segment(ode_id(nag,VACC_LOW,S),nag),E,I,densities.segment(ode_id(nag,VACC_LOW,R),nag)*/
        for(size_t i=0;i<nag;i++)
        {
            densities[ode_id(nag,LOW,E1,i)]=seeding_infectious[i]*Npop[i]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);
            densities[ode_id(nag,HIGH,E1,i)]=seeding_infectious[i]*Npop[i+nag]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);
            densities[ode_id(nag,PREG,E1,i)]=seeding_infectious[i]*Npop[i+2*nag]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);

            densities[ode_id(nag,LOW,S,i)]=Npop[i]-densities[ode_id(nag,LOW,E1,i)];
            densities[ode_id(nag,HIGH,S,i)]=Npop[i+nag]-densities[ode_id(nag,HIGH,E1,i)];
            densities[ode_id(nag,PREG,S,i)]=Npop[i+2*nag]-densities[ode_id(nag,PREG,E1,i)];
        }

        cases_t cases;

        cases.cases = Eigen::MatrixXd::Zero( (end_time-start_time).hours()
                /minimal_resolution, 
                contact_regular.cols()*group_types.size()/2);
        cases.times.reserve( (end_time-start_time).hours()
                /minimal_resolution );

        static bt::time_duration dt = bt::hours( 6 );
        while (cases.times.size()<cases.cases.rows())
        {
            auto next_time = current_time + bt::hours( minimal_resolution );

            Eigen::VectorXd vacc_rates;
            if (vaccine_programme.dates.size() > 0)
            {
                while (date_id < ((int)vaccine_programme.dates.size())-1 && 
                        current_time == vaccine_programme.dates[date_id+1] )
                {
                    ++date_id;
                }
            } else {
                // Legacy mode
                date_id=floor((current_time-start_time).hours()/24.0-44);
            }

            bool time_changed_for_vacc = false;
            if (date_id < ((int)vaccine_programme.dates.size())-1 && 
                        date_id >= -1 && 
                    next_time > vaccine_programme.dates[date_id+1] )
            {
                next_time = 
                    vaccine_programme.dates[date_id+1];
                time_changed_for_vacc = true;
            }

            if (date_id >= 0 &&
                    date_id < vaccine_programme.calendar.rows() )
                vacc_rates = vaccine_programme.calendar.row(date_id); 

            auto n_cases = new_cases( densities, current_time,
                    next_time, dt,
                    Npop,
                    vacc_rates,
                    vaccine_programme.efficacy_year,
                    transmission_regular,
                    a1, a2, g1, g2 );
            current_time = next_time;

            assert(step_count < cases.cases.rows());

            cases.cases.row(step_count) += n_cases;
            if (!time_changed_for_vacc) 
            {
                ++step_count;
                cases.times.push_back( current_time );
            }
        }
       
        return cases;
    } 

    Eigen::MatrixXd days_to_weeks_5AG(const cases_t &simulation)
    {

        size_t weeks =  (simulation.times.back() - simulation.times.front())
            .hours()/(24*7) + 1;
        auto result_days = simulation.cases;
        /*initialisation*/
        Eigen::MatrixXd result_weeks = 
            Eigen::MatrixXd::Zero( weeks, 5 );

        size_t j = 0;
        for(size_t i=0; i<weeks; i++)
        {
            auto startWeek = simulation.times[j];
            while( j < simulation.times.size() &&
                    (simulation.times[j]-startWeek).hours()/(24.0)<7.0 )
            {
                result_weeks(i,0)+=result_days(j,0)+result_days(j,1)+result_days(j,7)+result_days(j,8)+result_days(j,14)+result_days(j,15);
                result_weeks(i,1)+=result_days(j,2)+result_days(j,9)+result_days(j,16);
                result_weeks(i,2)+=result_days(j,3)+result_days(j,4)+result_days(j,10)+result_days(j,11)+result_days(j,17)+result_days(j,18);
                result_weeks(i,3)+=result_days(j,5)+result_days(j,12)+result_days(j,19);
                result_weeks(i,4)+=result_days(j,6)+result_days(j,13)+result_days(j,20);
                ++j;
            }
        }

        return result_weeks;
    }

    long double log_likelihood( double epsilon, double psi, 
            size_t predicted, double population_size, 
            int ili_cases, int ili_monitored,
            int confirmed_positive, int confirmed_samples, 
            int depth )
    {
        int Z_in_mon=(int)round(predicted*ili_monitored/population_size);
        int n=confirmed_samples;
        int m=ili_cases;

        /*h.init=max(n.plus-Z.in.mon,0)*/
        int h_init;
        if(confirmed_positive>Z_in_mon)
            h_init=confirmed_positive-Z_in_mon;
        else
            h_init=0;

        if(h_init>depth)
        {
            return -(h_init-depth)*
                std::numeric_limits<double>::max()/1e6;
        }

        /*define the first aij*/
        long double aij=pow(epsilon,confirmed_positive)*exp(-psi*ili_monitored*epsilon);

        if(confirmed_positive<n)
            for(int g=confirmed_positive; g<n; g++)
                aij*=m-g;

        if((Z_in_mon==confirmed_positive)&&(Z_in_mon>0))
            for(int g=1; g<=confirmed_positive;g++)
                aij*=g;

        if((confirmed_positive>0)&&(confirmed_positive<Z_in_mon))
            for(int g=0;g<confirmed_positive;g++)
                aij*=Z_in_mon-g;

        if((Z_in_mon>0)&&(confirmed_positive>Z_in_mon))
            for(int g=0; g<Z_in_mon;g++)
                aij*=confirmed_positive-g;

        if(confirmed_positive>Z_in_mon)
            aij*=pow(psi*ili_monitored,confirmed_positive-Z_in_mon);

        if(confirmed_positive<Z_in_mon)
            aij*=pow(1-epsilon,Z_in_mon-confirmed_positive);

        /*store the values of the first aij on the current line*/
        long double aij_seed=aij;
        int k_seed=confirmed_positive;

        auto likelihood_AG_week=aij;

        /*Calculation of the first line*/
        int max_m_plus;
        if(Z_in_mon+h_init>m-n+confirmed_positive)
            max_m_plus=m-n+confirmed_positive;
        else
            max_m_plus=Z_in_mon+h_init;

        if(max_m_plus>confirmed_positive)
            for(int k=confirmed_positive+1;k<=max_m_plus;k++)
            {
                aij*=(k*(m-k-n+confirmed_positive+1)*(Z_in_mon-k+h_init+1)*epsilon)/((m-k+1)*(k-confirmed_positive)*(k-h_init)*(1-epsilon));
                likelihood_AG_week+=aij;
            }

        /*top_sum=min(depth,m-n+confirmed_positive)*/
        int top_sum;
        if(depth<m-n+confirmed_positive)
            top_sum=depth;
        else
            top_sum=m-n+confirmed_positive;

        if(h_init<top_sum)
            for(int h=h_init+1;h<=top_sum;h++)
            {
                if(h>confirmed_positive) /*diagonal increment*/
                {
                    k_seed++;
                    aij_seed*=(k_seed*(m-k_seed-n+confirmed_positive+1)*psi*epsilon*ili_monitored)/((m-k_seed+1)*(k_seed-confirmed_positive)*h);
                }

                if(h<=confirmed_positive) /*vertical increment*/
                    aij_seed*=((k_seed-h+1)*psi*ili_monitored*(1-epsilon))/((Z_in_mon-k_seed+h)*h);

                aij=aij_seed;
                likelihood_AG_week+=aij;

                /*calculation of the line*/
                if(Z_in_mon+h>m-n+confirmed_positive)
                    max_m_plus=m-n+confirmed_positive;
                else
                    max_m_plus=Z_in_mon+h;

                if(max_m_plus>k_seed)
                    for(int k=k_seed+1;k<=max_m_plus;k++)
                    {
                        aij*=(k*(m-k-n+confirmed_positive+1)*(Z_in_mon-k+h+1)*epsilon)/((m-k+1)*(k-confirmed_positive)*(k-h)*(1-epsilon));
                        likelihood_AG_week+=aij;
                    }
            }

        auto ll = log(likelihood_AG_week);
        if (!std::isfinite(ll))
        {
            /*Rcpp::Rcerr << "Numerical error detected for week " 
              << week << " and age group " << i << std::endl;
              Rcpp::Rcerr << "Predicted number of cases is: "
              << result_by_week(week,i) << std::endl;*/
            ll = -std::numeric_limits<double>::max()/1e7;
        }
        return ll;
    }

    double log_likelihood_hyper_poisson(const std::vector<double> &eps, 
            double psi, const Eigen::MatrixXd &result_by_week,
            const Eigen::MatrixXi &ili, const Eigen::MatrixXi &mon_pop, 
            const Eigen::MatrixXi &n_pos, const Eigen::MatrixXi &n_samples, 
            //int * n_ILI, int * mon_popu, int * n_posi, int * n_sampled, 
            double * pop_5AG_RCGP, int depth)
    {
        int g, i, k, k_seed, week, h, h_init, top_sum;
        int Z_in_mon, n, m, n_plus, max_m_plus, pop_mon;

        long double result=0.0;
        for(int i=0;i<5;i++)
        {
            auto epsilon=eps[i];
            for(int week=0;week<result_by_week.rows();week++)
            {
                result += log_likelihood( epsilon, psi, 
                        result_by_week(week,i), pop_5AG_RCGP[i],
                        ili(week,i), mon_pop(week,i),
                        n_pos(week,i), n_samples(week,i), depth );
            }
            
        }
        return(result);
    }

    /// Return the log prior probability of the proposed parameters - current parameters
    //
    // \param susceptibility whether to use the prior based on 2003/04
    double log_prior( const parameter_set &proposed,
            const parameter_set &current,
            bool susceptibility ) {

        // Parameters should be valid
        if( 
                proposed.epsilon[0] <= 0 || proposed.epsilon[0] >= 1 ||
                proposed.epsilon[2] <= 0 || proposed.epsilon[2] >= 1 ||
                proposed.epsilon[4] <= 0 || proposed.epsilon[4] >= 1 ||
                proposed.psi < 0 || proposed.psi > 1 ||
                proposed.transmissibility < 0 ||
                proposed.susceptibility[0] < 0 || proposed.susceptibility[1] > 1 ||
                proposed.susceptibility[3] < 0 || proposed.susceptibility[3] > 1 ||
                proposed.susceptibility[6] < 0 || proposed.susceptibility[6] > 1 ||
                proposed.init_pop<log(0.00001) || proposed.init_pop>log(10)
          )
            return log(0);

        double log_prior = 0;
        if (!susceptibility)
        {
            /*Prior for the transmissibility; year other than 2003/04*/
            /*correction for a normal prior with mu=0.1653183 and sd=0.02773053*/
            /*prior on q*/
            log_prior=(current.transmissibility-proposed.transmissibility)*(current.transmissibility+proposed.transmissibility-0.3306366)*650.2099;
        } else {

            /*prior on the susceptibility (year 2003/04)*/

            /*correction for a normal prior with mu=0.688 and sd=0.083 for the 0-14 */
            log_prior=(current.susceptibility[0]-proposed.susceptibility[0])*(current.susceptibility[0]+proposed.susceptibility[0]-1.376)*145.1589/2;
            /*correction for a normal prior with mu=0.529 and sd=0.122 for the 15-64 */
            log_prior+=(current.susceptibility[3]-proposed.susceptibility[3])*(current.susceptibility[3]+proposed.susceptibility[3]-1.058)*67.18624/2;
            /*correction for a normal prior with mu=0.523 and sd=0.175 for the 65+ */
            log_prior+=(current.susceptibility[6]-proposed.susceptibility[6])*(current.susceptibility[6]+proposed.susceptibility[6]-1.046)*32.65306/2;
        }

        /*Prior for the ascertainment probabilities*/

        /*correct for the prior from serology season (lognormal):"0-14" lm=-4.493789, ls=0.2860455*/
        log_prior += log(current.epsilon[0])-log(proposed.epsilon[0])+(log(current.epsilon[0])-log(proposed.epsilon[0]))*(log(current.epsilon[0])+log(proposed.epsilon[0])+8.987578)*6.110824;

        /*correct for the prior from serology season (lognormal):"15-64" lm=-4.117028, ls=0.4751615*/
        log_prior += log(current.epsilon[2])-log(proposed.epsilon[2])+(log(current.epsilon[2])-log(proposed.epsilon[2]))*(log(current.epsilon[2])+log(proposed.epsilon[2])+8.234056)*2.21456;

        /*correct for the prior from serology season (lognormal):"65+" lm=-2.977965, ls=1.331832*/
        log_prior += log(current.epsilon[4])-log(proposed.epsilon[4])+(log(current.epsilon[4])-log(proposed.epsilon[4]))*(log(current.epsilon[4])+log(proposed.epsilon[4])+5.95593)*0.2818844;
        return log_prior;
    }

}
