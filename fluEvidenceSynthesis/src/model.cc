#include "model.hh"

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

    // Could convert this to template code if needed for performance
    struct seir_t 
    {
        // Note can also be used for delta
        Eigen::VectorXd s;
        Eigen::VectorXd e1;
        Eigen::VectorXd e2;
        Eigen::VectorXd i1;
        Eigen::VectorXd i2;
        Eigen::VectorXd r;

        seir_t() {}
        seir_t( size_t dim )
        {
            /*s = Eigen::VectorXd( dim );
            e1 = Eigen::VectorXd( dim );
            e2 = Eigen::VectorXd( dim );
            i1 = Eigen::VectorXd( dim );
            i2 = Eigen::VectorXd( dim );
            r = Eigen::VectorXd( dim );*/
            s = Eigen::VectorXd::Zero( dim );
            e1 = Eigen::VectorXd::Zero( dim );
            e2 = Eigen::VectorXd::Zero( dim );
            i1 = Eigen::VectorXd::Zero( dim );
            i2 = Eigen::VectorXd::Zero( dim );
            r = Eigen::VectorXd::Zero( dim );
        }
    };

    enum group_type_t { LOW, HIGH, PREG,
        VACC_LOW, VACC_HIGH, VACC_PREG };
    std::vector<group_type_t> group_types = { LOW, HIGH, PREG,
        VACC_LOW, VACC_HIGH, VACC_PREG };

    inline Eigen::VectorXd new_cases( 
            std::vector<seir_t> &densities,
            const boost::posix_time::ptime &start_time,
            const boost::posix_time::ptime &end_time, 
            const boost::posix_time::time_duration &dt,
            const std::vector<double> &Npop,
            const Eigen::MatrixXd &vaccine_rates, // If empty, rate of zero is assumed
            const std::array<double,7> &vaccine_efficacy_year,
            const Eigen::MatrixXd &transmission_regular,
            double a1, double a2, double g1, double g2
            )
    {
        namespace bt = boost::posix_time;

        double h_step = dt.hours()/24.0;

        auto nag = transmission_regular.cols();
        Eigen::VectorXd results = Eigen::VectorXd::Zero(nag*3);

        std::vector<seir_t> deltas;
        for( auto &gt : group_types ) {
            deltas.push_back( seir_t(nag) );
        }

        for(bt::ptime current_time = start_time; 
                current_time < end_time;)
        {
            current_time += dt;
            if ((current_time > end_time))
            {
                h_step = (current_time - end_time).hours()/24.0;
                current_time = end_time;
            }

            for(size_t i=0;i<nag;i++)
            {
                /*rate of depletion of susceptible*/
                deltas[VACC_LOW].s[i]=0;
                for(size_t j=0;j<nag;j++)
                    deltas[VACC_LOW].s[i]+=transmission_regular(i,j)*(densities[VACC_LOW].i1[j]+densities[VACC_LOW].i2[j]+densities[VACC_HIGH].i1[j]+densities[VACC_HIGH].i2[j]+densities[VACC_PREG].i1[j]+densities[VACC_PREG].i2[j]+densities[LOW].i1[j]+densities[LOW].i2[j]+densities[HIGH].i1[j]+densities[HIGH].i2[j]+densities[PREG].i1[j]+densities[PREG].i2[j]);

                deltas[VACC_HIGH].s[i]=deltas[VACC_LOW].s[i];
                deltas[VACC_PREG].s[i]=deltas[VACC_LOW].s[i];
                deltas[LOW].s[i]=deltas[VACC_LOW].s[i];
                deltas[HIGH].s[i]=deltas[VACC_LOW].s[i];
                deltas[PREG].s[i]=deltas[VACC_LOW].s[i];

                deltas[VACC_LOW].s[i]*=-densities[VACC_LOW].s[i];
                deltas[VACC_HIGH].s[i]*=-densities[VACC_HIGH].s[i];
                deltas[VACC_PREG].s[i]*=-densities[VACC_PREG].s[i];
                deltas[LOW].s[i]*=-densities[LOW].s[i];
                deltas[HIGH].s[i]*=-densities[HIGH].s[i];
                deltas[PREG].s[i]*=-densities[PREG].s[i];
            }

            /*rate of passing between states of infection*/
            for ( auto &gt : group_types)
            {
                deltas[gt].e1=-deltas[gt].s-a1*densities[gt].e1;
                deltas[gt].e2=a1*densities[gt].e1-a2*densities[gt].e2;

                deltas[gt].i1=a2*densities[gt].e2-g1*densities[gt].i1;
                deltas[gt].i2=g1*densities[gt].i1-g2*densities[gt].i2;
                deltas[gt].r=g2*densities[gt].i2;
            }

            /*Vaccine bit*/

            if ( vaccine_rates.size() > 0 )
            {
                for(size_t i=0;i<nag;i++)
                {
                    double vacc_prov=Npop[i]*vaccine_rates(i)/(densities[LOW].s[i]+densities[LOW].e1[i]+densities[LOW].e2[i]+densities[LOW].i1[i]+densities[LOW].i2[i]+densities[LOW].r[i]);
                    /*surv[i]+=vaccination_calendar[cal_time*21+i];*/
                    double vacc_prov_r=Npop[i+nag]*vaccine_rates(i+nag)/(densities[HIGH].s[i]+densities[HIGH].e1[i]+densities[HIGH].e2[i]+densities[HIGH].i1[i]+densities[HIGH].i2[i]+densities[HIGH].r[i]);
                    double vacc_prov_p=0; /*Npop[i+2*nag]*vaccination_calendar[cal_time*21+i+2*nag]/(densities[PREG].s[i]+densities[PREG].e1[i]+densities[PREG].e2[i]+densities[PREG].i1[i]+densities[PREG].i2[i]+densities[PREG].r[i]);*/

                    deltas[VACC_LOW].s[i]+=densities[LOW].s[i]*vacc_prov*(1-vaccine_efficacy_year[i]);
                    deltas[VACC_HIGH].s[i]+=densities[HIGH].s[i]*vacc_prov_r*(1-vaccine_efficacy_year[i]);
                    deltas[VACC_PREG].s[i]+=densities[PREG].s[i]*vacc_prov_p*(1-vaccine_efficacy_year[i]);
                    deltas[LOW].s[i]-=densities[LOW].s[i]*vacc_prov;
                    deltas[HIGH].s[i]-=densities[HIGH].s[i]*vacc_prov_r;
                    deltas[PREG].s[i]-=densities[PREG].s[i]*vacc_prov_p;

                    deltas[VACC_LOW].e1[i]+=densities[LOW].e1[i]*vacc_prov;
                    deltas[VACC_HIGH].e1[i]+=densities[HIGH].e1[i]*vacc_prov_r;
                    deltas[VACC_PREG].e1[i]+=densities[PREG].e1[i]*vacc_prov_p;
                    deltas[LOW].e1[i]-=densities[LOW].e1[i]*vacc_prov;
                    deltas[HIGH].e1[i]-=densities[HIGH].e1[i]*vacc_prov_r;
                    deltas[PREG].e1[i]-=densities[PREG].e1[i]*vacc_prov_p;

                    deltas[VACC_LOW].e2[i]+=densities[LOW].e2[i]*vacc_prov;
                    deltas[VACC_HIGH].e2[i]+=densities[HIGH].e2[i]*vacc_prov_r;
                    deltas[VACC_PREG].e2[i]+=densities[PREG].e2[i]*vacc_prov_p;
                    deltas[LOW].e2[i]-=densities[LOW].e2[i]*vacc_prov;
                    deltas[HIGH].e2[i]-=densities[HIGH].e2[i]*vacc_prov_r;
                    deltas[PREG].e2[i]-=densities[PREG].e2[i]*vacc_prov_p;

                    deltas[VACC_LOW].i1[i]+=densities[LOW].i1[i]*vacc_prov;
                    deltas[VACC_HIGH].i1[i]+=densities[HIGH].i1[i]*vacc_prov_r;
                    deltas[VACC_PREG].i1[i]+=densities[PREG].i1[i]*vacc_prov_p;
                    deltas[LOW].i1[i]-=densities[LOW].i1[i]*vacc_prov;
                    deltas[HIGH].i1[i]-=densities[HIGH].i1[i]*vacc_prov_r;
                    deltas[PREG].i1[i]-=densities[PREG].i1[i]*vacc_prov_p;

                    deltas[VACC_LOW].i2[i]+=densities[LOW].i2[i]*vacc_prov;
                    deltas[VACC_HIGH].i2[i]+=densities[HIGH].i2[i]*vacc_prov_r;
                    deltas[VACC_PREG].i2[i]+=densities[PREG].i2[i]*vacc_prov_p;
                    deltas[LOW].i2[i]-=densities[LOW].i2[i]*vacc_prov;
                    deltas[HIGH].i2[i]-=densities[HIGH].i2[i]*vacc_prov_r;
                    deltas[PREG].i2[i]-=densities[PREG].i2[i]*vacc_prov_p;

                    deltas[VACC_LOW].r[i]+=densities[LOW].r[i]*vacc_prov+densities[LOW].s[i]*vacc_prov*vaccine_efficacy_year[i];
                    deltas[VACC_HIGH].r[i]+=densities[HIGH].r[i]*vacc_prov_r+densities[HIGH].s[i]*vacc_prov_r*vaccine_efficacy_year[i];
                    deltas[VACC_PREG].r[i]+=densities[PREG].r[i]*vacc_prov_p+densities[PREG].s[i]*vacc_prov_p*vaccine_efficacy_year[i];
                    deltas[LOW].r[i]-=densities[LOW].r[i]*vacc_prov;
                    deltas[HIGH].r[i]-=densities[HIGH].r[i]*vacc_prov_r;
                    deltas[PREG].r[i]-=densities[PREG].r[i]*vacc_prov_p;
                }
            }

            /*update the different classes*/
            for( auto &gt : group_types ) {
                densities[gt].s += h_step * deltas[gt].s; 
                densities[gt].e1 += h_step * deltas[gt].e1; 
                densities[gt].e2 += h_step * deltas[gt].e2; 
                densities[gt].i1 += h_step * deltas[gt].i1; 
                densities[gt].i2 += h_step * deltas[gt].i2; 
                densities[gt].r += h_step * deltas[gt].r; 
            }

            results.block( 0, 0, nag, 1 ) += a2*(densities[VACC_LOW].e2+densities[LOW].e2)*h_step;
            results.block( nag, 0, nag, 1 ) += a2*(densities[VACC_HIGH].e2+densities[HIGH].e2)*h_step;
            results.block( 2*nag, 0, nag, 1 ) += a2*(densities[VACC_PREG].e2+densities[PREG].e2)*h_step;
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

        std::vector<seir_t> densities;
        std::vector<seir_t> deltas;
        for( auto &gt : group_types ) {
            densities.push_back( seir_t(nag) );
            deltas.push_back( seir_t(nag) );
        }

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
        for(size_t i=0;i<transmission_regular.rows();i++)
        {
            for(size_t j=0;j<transmission_regular.cols();j++) {
                transmission_regular(i,j)*=q*s_profile[i];
            }
        }

        /*initialisation, densities[VACC_LOW].s,E,I,densities[VACC_LOW].r*/
        for(size_t i=0;i<nag;i++)
        {
            densities[LOW].e1[i]=seeding_infectious[i]*Npop[i]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);
            densities[HIGH].e1[i]=seeding_infectious[i]*Npop[i+nag]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);
            densities[PREG].e1[i]=seeding_infectious[i]*Npop[i+2*nag]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);

            densities[LOW].s[i]=Npop[i]-densities[LOW].e1[i];
            densities[HIGH].s[i]=Npop[i+nag]-densities[HIGH].e1[i];
            densities[PREG].s[i]=Npop[i+2*nag]-densities[PREG].e1[i];
        }

        cases_t cases;

        cases.cases = Eigen::MatrixXd::Zero( (end_time-start_time).hours()
                /minimal_resolution, 
                contact_regular.cols()*group_types.size()/2);
        cases.times.reserve( (end_time-start_time).hours()
                /minimal_resolution );

        /*Eigen::VectorXd new_cases( 
            std::vector<seir_t> &densities,
            const boost::posix_time::ptime &start_time,
            const boost::posix_time::ptime &end_time, 
            const boost::posix_time::time_duration &dt,
            const std::vector<double> &Npop,
            const Eigen::MatrixXd &vaccine_rates, // If empty, rate of zero is assumed
            const Eigen::VectorXd &vaccine_efficacy_year,
            const Eigen::MatrixXd &transmission_regular,
            double a1, double a2, double g1, double g2
            )*/

        while (current_time < end_time)
        {
            bt::time_duration dt = bt::hours( 6 );
            auto next_time = current_time + bt::hours( minimal_resolution );

            bool time_changed_for_vacc = false;
            if (date_id < ((int)vaccine_programme.dates.size())-1 && 
                        date_id >= -1 && 
                    next_time > vaccine_programme.dates[date_id+1] )
            {
                next_time = 
                    vaccine_programme.dates[date_id+1];
                dt = next_time - current_time;
                time_changed_for_vacc = true;
            }

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
    cases_t one_year_SEIR_with_vaccination2(
            const std::vector<double> &Npop, double * seeding_infectious, 
            const double tlatent, const double tinfectious, 
            const std::vector<double> &s_profile, 
            const Eigen::MatrixXd &contact_regular, double q,
            const vaccine::vaccine_t &vaccine_programme,
            size_t minimal_resolution )
    {
        namespace bt = boost::posix_time;

        const size_t nag = contact_regular.rows(); // No. of age groups

        std::vector<seir_t> densities;
        std::vector<seir_t> deltas;
        for( auto &gt : group_types ) {
            densities.push_back( seir_t(nag) );
            deltas.push_back( seir_t(nag) );
        }

        double vacc_prov, vacc_prov_p, vacc_prov_r;
        int i, j;
        double a1, a2, g1, g2 /*, surv[7]={0,0,0,0,0,0,0}*/;

        double h_step = 0.25;
        bt::time_duration dt = bt::hours( 24 * h_step );

        // We start at week 35. Week 1 is the first week that ends in this year
        auto current_time = getTimeFromWeekYear( 35, 1970 );
        if (vaccine_programme.dates.size()!=0)
            current_time = getTimeFromWeekYear( 35, 
                    vaccine_programme.dates[0].date().year() );

        auto start_time = current_time;
        auto end_time = current_time + bt::hours(364*24);
        auto last_recorded = current_time;
        size_t step_count = 0;

        a1=2/tlatent;
        a2=a1;
        g1=2/tinfectious;
        g2=g1;

        int date_id = -1;

        /*initialisation, transmission matrix*/
        Eigen::MatrixXd transmission_regular(contact_regular);
        for(i=0;i<transmission_regular.rows();i++)
        {
            for(size_t j=0;j<transmission_regular.cols();j++) {
                transmission_regular(i,j)*=q*s_profile[i];
            }
        }

        Eigen::VectorXd total_of_new_cases_per_day = 
            Eigen::VectorXd::Zero( nag );
        Eigen::VectorXd total_of_new_cases_per_day_r = 
            Eigen::VectorXd::Zero( nag );
        Eigen::VectorXd total_of_new_cases_per_day_p = 
            Eigen::VectorXd::Zero( nag );

        /*initialisation, densities[VACC_LOW].s,E,I,densities[VACC_LOW].r*/
        for(i=0;i<nag;i++)
        {
            densities[LOW].e1[i]=seeding_infectious[i]*Npop[i]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);
            densities[HIGH].e1[i]=seeding_infectious[i]*Npop[i+nag]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);
            densities[PREG].e1[i]=seeding_infectious[i]*Npop[i+2*nag]/(Npop[i]+Npop[i+nag]+Npop[i+2*nag]);

            densities[LOW].s[i]=Npop[i]-densities[LOW].e1[i];
            densities[HIGH].s[i]=Npop[i+nag]-densities[HIGH].e1[i];
            densities[PREG].s[i]=Npop[i+2*nag]-densities[PREG].e1[i];
        }

        cases_t cases;

        cases.cases = Eigen::MatrixXd( (end_time-start_time).hours()
                /minimal_resolution, 
                contact_regular.cols()*group_types.size()/2);
        cases.times.resize( (end_time-start_time).hours()
                /minimal_resolution );

        for(; current_time <= end_time; current_time += dt)
        {
            // Check we don't overstep, either vaccine calendar next date
            // or minimal_resolution
            if ((current_time - last_recorded).hours()>minimal_resolution)
            {
                auto new_current_time = 
                    last_recorded + bt::hours(minimal_resolution);
                h_step = (dt-(current_time-new_current_time)).hours()/24.0;
                current_time = new_current_time;
            }
            if (date_id < ((int)vaccine_programme.dates.size())-1 && 
                        date_id >= -1 && 
                    current_time > vaccine_programme.dates[date_id+1] )
            {
                auto new_current_time = 
                    vaccine_programme.dates[date_id+1];
                h_step = (dt-(current_time-new_current_time)).hours()/24.0;
                current_time = new_current_time;
            } else {
                h_step = dt.hours()/24.0;
            }

            if((current_time - last_recorded).hours()==minimal_resolution)
            {
                for(i=0;i<nag;i++)
                {
                    cases.cases(step_count,i)=total_of_new_cases_per_day[i];
                    cases.cases(step_count,nag+i)=total_of_new_cases_per_day_r[i];
                    cases.cases(step_count,2*nag+i)=total_of_new_cases_per_day_p[i];
                }
                cases.times[step_count] = current_time;
                last_recorded = current_time;
                //cases.times[(int)t] = dt;

                total_of_new_cases_per_day = Eigen::VectorXd::Zero( nag );
                total_of_new_cases_per_day_r = Eigen::VectorXd::Zero( nag );
                total_of_new_cases_per_day_p = Eigen::VectorXd::Zero( nag );
                ++step_count;
            }
            
            for(i=0;i<nag;i++)
            {
                /*rate of depletion of susceptible*/
                deltas[VACC_LOW].s[i]=0;
                for(j=0;j<nag;j++)
                    deltas[VACC_LOW].s[i]+=transmission_regular(i,j)*(densities[VACC_LOW].i1[j]+densities[VACC_LOW].i2[j]+densities[VACC_HIGH].i1[j]+densities[VACC_HIGH].i2[j]+densities[VACC_PREG].i1[j]+densities[VACC_PREG].i2[j]+densities[LOW].i1[j]+densities[LOW].i2[j]+densities[HIGH].i1[j]+densities[HIGH].i2[j]+densities[PREG].i1[j]+densities[PREG].i2[j]);

                deltas[VACC_HIGH].s[i]=deltas[VACC_LOW].s[i];
                deltas[VACC_PREG].s[i]=deltas[VACC_LOW].s[i];
                deltas[LOW].s[i]=deltas[VACC_LOW].s[i];
                deltas[HIGH].s[i]=deltas[VACC_LOW].s[i];
                deltas[PREG].s[i]=deltas[VACC_LOW].s[i];

                deltas[VACC_LOW].s[i]*=-densities[VACC_LOW].s[i];
                deltas[VACC_HIGH].s[i]*=-densities[VACC_HIGH].s[i];
                deltas[VACC_PREG].s[i]*=-densities[VACC_PREG].s[i];
                deltas[LOW].s[i]*=-densities[LOW].s[i];
                deltas[HIGH].s[i]*=-densities[HIGH].s[i];
                deltas[PREG].s[i]*=-densities[PREG].s[i];
            }

            /*rate of passing between states of infection*/
            for ( auto &gt : group_types)
            {
                deltas[gt].e1=-deltas[gt].s-a1*densities[gt].e1;
                deltas[gt].e2=a1*densities[gt].e1-a2*densities[gt].e2;

                deltas[gt].i1=a2*densities[gt].e2-g1*densities[gt].i1;
                deltas[gt].i2=g1*densities[gt].i1-g2*densities[gt].i2;
                deltas[gt].r=g2*densities[gt].i2;
            }

            /*Vaccine bit*/

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

            if (date_id >= 0 &&
                    date_id < vaccine_programme.calendar.rows() )
            {
                for(i=0;i<nag;i++)
                {
                    vacc_prov=Npop[i]*vaccine_programme.calendar(date_id,i)/(densities[LOW].s[i]+densities[LOW].e1[i]+densities[LOW].e2[i]+densities[LOW].i1[i]+densities[LOW].i2[i]+densities[LOW].r[i]);
                    /*surv[i]+=vaccination_calendar[cal_time*21+i];*/
                    vacc_prov_r=Npop[i+nag]*vaccine_programme.calendar(date_id,i+nag)/(densities[HIGH].s[i]+densities[HIGH].e1[i]+densities[HIGH].e2[i]+densities[HIGH].i1[i]+densities[HIGH].i2[i]+densities[HIGH].r[i]);
                    vacc_prov_p=0; /*Npop[i+2*nag]*vaccination_calendar[cal_time*21+i+2*nag]/(densities[PREG].s[i]+densities[PREG].e1[i]+densities[PREG].e2[i]+densities[PREG].i1[i]+densities[PREG].i2[i]+densities[PREG].r[i]);*/

                    deltas[VACC_LOW].s[i]+=densities[LOW].s[i]*vacc_prov*(1-vaccine_programme.efficacy_year[i]);
                    deltas[VACC_HIGH].s[i]+=densities[HIGH].s[i]*vacc_prov_r*(1-vaccine_programme.efficacy_year[i]);
                    deltas[VACC_PREG].s[i]+=densities[PREG].s[i]*vacc_prov_p*(1-vaccine_programme.efficacy_year[i]);
                    deltas[LOW].s[i]-=densities[LOW].s[i]*vacc_prov;
                    deltas[HIGH].s[i]-=densities[HIGH].s[i]*vacc_prov_r;
                    deltas[PREG].s[i]-=densities[PREG].s[i]*vacc_prov_p;

                    deltas[VACC_LOW].e1[i]+=densities[LOW].e1[i]*vacc_prov;
                    deltas[VACC_HIGH].e1[i]+=densities[HIGH].e1[i]*vacc_prov_r;
                    deltas[VACC_PREG].e1[i]+=densities[PREG].e1[i]*vacc_prov_p;
                    deltas[LOW].e1[i]-=densities[LOW].e1[i]*vacc_prov;
                    deltas[HIGH].e1[i]-=densities[HIGH].e1[i]*vacc_prov_r;
                    deltas[PREG].e1[i]-=densities[PREG].e1[i]*vacc_prov_p;

                    deltas[VACC_LOW].e2[i]+=densities[LOW].e2[i]*vacc_prov;
                    deltas[VACC_HIGH].e2[i]+=densities[HIGH].e2[i]*vacc_prov_r;
                    deltas[VACC_PREG].e2[i]+=densities[PREG].e2[i]*vacc_prov_p;
                    deltas[LOW].e2[i]-=densities[LOW].e2[i]*vacc_prov;
                    deltas[HIGH].e2[i]-=densities[HIGH].e2[i]*vacc_prov_r;
                    deltas[PREG].e2[i]-=densities[PREG].e2[i]*vacc_prov_p;

                    deltas[VACC_LOW].i1[i]+=densities[LOW].i1[i]*vacc_prov;
                    deltas[VACC_HIGH].i1[i]+=densities[HIGH].i1[i]*vacc_prov_r;
                    deltas[VACC_PREG].i1[i]+=densities[PREG].i1[i]*vacc_prov_p;
                    deltas[LOW].i1[i]-=densities[LOW].i1[i]*vacc_prov;
                    deltas[HIGH].i1[i]-=densities[HIGH].i1[i]*vacc_prov_r;
                    deltas[PREG].i1[i]-=densities[PREG].i1[i]*vacc_prov_p;

                    deltas[VACC_LOW].i2[i]+=densities[LOW].i2[i]*vacc_prov;
                    deltas[VACC_HIGH].i2[i]+=densities[HIGH].i2[i]*vacc_prov_r;
                    deltas[VACC_PREG].i2[i]+=densities[PREG].i2[i]*vacc_prov_p;
                    deltas[LOW].i2[i]-=densities[LOW].i2[i]*vacc_prov;
                    deltas[HIGH].i2[i]-=densities[HIGH].i2[i]*vacc_prov_r;
                    deltas[PREG].i2[i]-=densities[PREG].i2[i]*vacc_prov_p;

                    deltas[VACC_LOW].r[i]+=densities[LOW].r[i]*vacc_prov+densities[LOW].s[i]*vacc_prov*vaccine_programme.efficacy_year[i];
                    deltas[VACC_HIGH].r[i]+=densities[HIGH].r[i]*vacc_prov_r+densities[HIGH].s[i]*vacc_prov_r*vaccine_programme.efficacy_year[i];
                    deltas[VACC_PREG].r[i]+=densities[PREG].r[i]*vacc_prov_p+densities[PREG].s[i]*vacc_prov_p*vaccine_programme.efficacy_year[i];
                    deltas[LOW].r[i]-=densities[LOW].r[i]*vacc_prov;
                    deltas[HIGH].r[i]-=densities[HIGH].r[i]*vacc_prov_r;
                    deltas[PREG].r[i]-=densities[PREG].r[i]*vacc_prov_p;
                }
            }

            /*update the different classes*/
            for( auto &gt : group_types ) {
                densities[gt].s += h_step * deltas[gt].s; 
                densities[gt].e1 += h_step * deltas[gt].e1; 
                densities[gt].e2 += h_step * deltas[gt].e2; 
                densities[gt].i1 += h_step * deltas[gt].i1; 
                densities[gt].i2 += h_step * deltas[gt].i2; 
                densities[gt].r += h_step * deltas[gt].r; 
            }

            total_of_new_cases_per_day += a2*(densities[VACC_LOW].e2+densities[LOW].e2)*h_step;
            total_of_new_cases_per_day_r += a2*(densities[VACC_HIGH].e2+densities[HIGH].e2)*h_step;
            total_of_new_cases_per_day_p += a2*(densities[VACC_PREG].e2+densities[PREG].e2)*h_step;
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

    double log_likelihood_hyper_poisson(const std::vector<double> &eps, 
            double psi, const Eigen::MatrixXd &result_by_week,
            const Eigen::MatrixXi &ili, const Eigen::MatrixXi &mon_pop, 
            const Eigen::MatrixXi &n_pos, const Eigen::MatrixXi &n_samples, 
            //int * n_ILI, int * mon_popu, int * n_posi, int * n_sampled, 
            double * pop_5AG_RCGP, int depth)
    {
        int g, i, k, k_seed, week, h, h_init, top_sum;
        int Z_in_mon, n, m, n_plus, max_m_plus, pop_mon;
        double result, epsilon;
        long double aij, aij_seed, likelihood_AG_week;

        result=0.0;
        for(i=0;i<5;i++)
        {
            epsilon=eps[i];
            for(week=0;week<result_by_week.rows();week++)
            {
                pop_mon=mon_pop(week,i);
                Z_in_mon=(int)round(result_by_week(week,i)*pop_mon/pop_5AG_RCGP[i]);
                n=n_samples(week,i);
                n_plus=n_pos(week,i);
                m=ili(week,i);

                /*h.init=max(n.plus-Z.in.mon,0)*/
                if(n_plus>Z_in_mon)
                    h_init=n_plus-Z_in_mon;
                else
                    h_init=0;

                if(h_init>depth)
                    return(-10000);

                /*define the first aij*/
                aij=pow(epsilon,n_plus)*exp(-psi*pop_mon*epsilon);

                if(n_plus<n)
                    for(g=n_plus; g<n; g++)
                        aij*=m-g;

                if((Z_in_mon==n_plus)&&(Z_in_mon>0))
                    for(g=1; g<=n_plus;g++)
                        aij*=g;

                if((n_plus>0)&&(n_plus<Z_in_mon))
                    for(g=0;g<n_plus;g++)
                        aij*=Z_in_mon-g;

                if((Z_in_mon>0)&&(n_plus>Z_in_mon))
                    for(g=0; g<Z_in_mon;g++)
                        aij*=n_plus-g;

                if(n_plus>Z_in_mon)
                    aij*=pow(psi*pop_mon,n_plus-Z_in_mon);

                if(n_plus<Z_in_mon)
                    aij*=pow(1-epsilon,Z_in_mon-n_plus);

                /*store the values of the first aij on the current line*/
                aij_seed=aij;
                k_seed=n_plus;

                likelihood_AG_week=aij;

                /*Calculation of the first line*/
                if(Z_in_mon+h_init>m-n+n_plus)
                    max_m_plus=m-n+n_plus;
                else
                    max_m_plus=Z_in_mon+h_init;

                if(max_m_plus>n_plus)
                    for(k=n_plus+1;k<=max_m_plus;k++)
                    {
                        aij*=(k*(m-k-n+n_plus+1)*(Z_in_mon-k+h_init+1)*epsilon)/((m-k+1)*(k-n_plus)*(k-h_init)*(1-epsilon));
                        likelihood_AG_week+=aij;
                    }

                /*top_sum=min(depth,m-n+n_plus)*/
                if(depth<m-n+n_plus)
                    top_sum=depth;
                else
                    top_sum=m-n+n_plus;

                if(h_init<top_sum)
                    for(h=h_init+1;h<=top_sum;h++)
                    {
                        if(h>n_plus) /*diagonal increment*/
                        {
                            k_seed++;
                            aij_seed*=(k_seed*(m-k_seed-n+n_plus+1)*psi*epsilon*pop_mon)/((m-k_seed+1)*(k_seed-n_plus)*h);
                        }

                        if(h<=n_plus) /*vertical increment*/
                            aij_seed*=((k_seed-h+1)*psi*pop_mon*(1-epsilon))/((Z_in_mon-k_seed+h)*h);

                        aij=aij_seed;
                        likelihood_AG_week+=aij;

                        /*calculation of the line*/
                        if(Z_in_mon+h>m-n+n_plus)
                            max_m_plus=m-n+n_plus;
                        else
                            max_m_plus=Z_in_mon+h;

                        if(max_m_plus>k_seed)
                            for(k=k_seed+1;k<=max_m_plus;k++)
                            {
                                aij*=(k*(m-k-n+n_plus+1)*(Z_in_mon-k+h+1)*epsilon)/((m-k+1)*(k-n_plus)*(k-h)*(1-epsilon));
                                likelihood_AG_week+=aij;
                            }
                    }

                result+=log(likelihood_AG_week);
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

};
