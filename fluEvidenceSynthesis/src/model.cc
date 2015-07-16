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

        void setToZero()
        {
            for (size_t i = 0; i < s.size(); ++i) {
                s[i] = 0;
                e1[i] = 0;
                e2[i] = 0;
                i1[i] = 0;
                i2[i] = 0;
                r[i] = 0;
            }
            /*s = Eigen::VectorXd::Constant(s.size(), 0);
            e1 = Eigen::VectorXd::Zero( e1.size() );
            e2 = Eigen::VectorXd::Zero( e2.size() );
            i1 = Eigen::VectorXd::Zero( i1.size() );
            i2 = Eigen::VectorXd::Zero( i2.size() );
            r = Eigen::VectorXd::Zero( r.size() );*/
         }
    };

    enum group_type_t { LOW, HIGH, PREG,
        VACC_LOW, VACC_HIGH, VACC_PREG };
    std::vector<group_type_t> group_types = { LOW, HIGH, PREG,
        VACC_LOW, VACC_HIGH, VACC_PREG };
 
    void one_year_SEIR_with_vaccination(double * result, const std::vector<double> &Npop, double * seeding_infectious, const double tlatent, const double tinfectious, const std::vector<double> &s_profile, const Eigen::MatrixXd &contact_regular, double q,
            const vaccine::vaccine_t &vaccine_programme )
    {
        namespace bt = boost::posix_time;

        const size_t nag = contact_regular.rows(); // No. of age groups

        static std::vector<seir_t> densities;
        static std::vector<seir_t> deltas;
        if (densities.size() == 0)
        {
            for( auto &gt : group_types ) {
                densities.push_back( seir_t(nag) );
                deltas.push_back( seir_t(nag) );
            }
        } else {
            for( auto &gt : group_types ) {
                densities[gt].setToZero();
                deltas[gt].setToZero();
            }
        }

        double vacc_prov, vacc_prov_p, vacc_prov_r;
        int i, j;
        double a1, a2, g1, g2, t /*, surv[7]={0,0,0,0,0,0,0}*/;
        int step_rate;

        step_rate=(int)(1/h_step);
        bt::time_duration dt = bt::hours( 24 * h_step );

        // We start at week 35. Week 1 is the first week that ends in this year
        auto current_time = getTimeFromWeekYear( 35, 1970 );
        if (vaccine_programme.dates.size()!=0)
            current_time = getTimeFromWeekYear( 35, 
                    vaccine_programme.dates[0].date().year() );

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

        for(t=0; t<no_days; t+=h_step)
        {
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
                        current_time >= vaccine_programme.dates[date_id+1] )
                {
                    ++date_id;
                }
            } else {
                // Legacy mode
                date_id=(int)(t)-44;
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

            if((((int)(t*step_rate))%step_rate)==step_rate/2)
            {
                for(i=0;i<nag;i++)
                {
                    result[(int)t*3*nag+i]=total_of_new_cases_per_day[i];
                    result[(int)t*3*nag+nag+i]=total_of_new_cases_per_day_r[i];
                    result[(int)t*3*nag+2*nag+i]=total_of_new_cases_per_day_p[i];
                }

                total_of_new_cases_per_day = Eigen::VectorXd::Zero( nag );
                total_of_new_cases_per_day_r = Eigen::VectorXd::Zero( nag );
                total_of_new_cases_per_day_p = Eigen::VectorXd::Zero( nag );
            }
            current_time += dt;
        }
    }

    void one_year_SEIR_without_vaccination(double * result, const std::vector<double> &Npop, double * seeding_infectious, const double tlatent, const double tinfectious, const std::vector<double> & s_profile, const Eigen::MatrixXd &contact_regular, double q)
    {
        double transmission_regular[NAG2];
        double S[7], E1[7], E2[7], I1[7], I2[7], R[7];
        double Sr[7], E1r[7], E2r[7], I1r[7], I2r[7], Rr[7];
        double Sp[7], E1p[7], E2p[7], I1p[7], I2p[7], Rp[7];
        double deltaS[7], deltaE1[7], deltaE2[7], deltaI1[7], deltaI2[7], deltaR[7];
        double deltaSr[7], deltaE1r[7], deltaE2r[7], deltaI1r[7], deltaI2r[7], deltaRr[7];
        double deltaSp[7], deltaE1p[7], deltaE2p[7], deltaI1p[7], deltaI2p[7], deltaRp[7];
        double total_of_new_cases_per_day[7], total_of_new_cases_per_day_r[7], total_of_new_cases_per_day_p[7];
        int i, j;
        double a1, a2, g1, g2, t;
        int step_rate;

        step_rate=(int)(1/h_step);

        a1=2/tlatent;
        a2=a1;
        g1=2/tinfectious;
        g2=g1;

        /*initialisation, transmission matrix*/
        /*initialisation, transmission matrix*/
        for(i=0;i<NAG;i++)
        {
            for(size_t j=0;j<NAG;j++) {
                transmission_regular[i*NAG+j]=q*contact_regular(i,j)*s_profile[i];
            }
        }

        /*initialisation, S,E,I,R*/
        for(i=0;i<NAG;i++)
        {
            E1[i]=seeding_infectious[i]*Npop[i]/(Npop[i]+Npop[i+NAG]+Npop[i+2*NAG]);
            E1r[i]=seeding_infectious[i]*Npop[i+NAG]/(Npop[i]+Npop[i+NAG]+Npop[i+2*NAG]);
            E1p[i]=seeding_infectious[i]*Npop[i+2*NAG]/(Npop[i]+Npop[i+NAG]+Npop[i+2*NAG]);
            E2[i]=0;
            E2r[i]=0;
            E2p[i]=0;
            I1[i]=0;
            I1r[i]=0;
            I1p[i]=0;
            I2[i]=0;
            I2r[i]=0;
            I2p[i]=0;
            R[i]=0;
            Rr[i]=0;
            Rp[i]=0;
            S[i]=Npop[i]-E1[i];
            Sr[i]=Npop[i+NAG]-E1r[i];
            Sp[i]=Npop[i+2*NAG]-E1p[i];

            total_of_new_cases_per_day[i]=0;
            total_of_new_cases_per_day_r[i]=0;
            total_of_new_cases_per_day_p[i]=0;
        }

        for(t=0; t<no_days; t+=h_step)
        {
            for(i=0;i<NAG;i++)
            {
                /*initialisation of the different increment for the euler algorithm*/
                deltaS[i]=0;
                deltaE1[i]=0;
                deltaE2[i]=0;
                deltaI1[i]=0;
                deltaI2[i]=0;
                deltaR[i]=0;

                deltaSr[i]=0;
                deltaE1r[i]=0;
                deltaE2r[i]=0;
                deltaI1r[i]=0;
                deltaI2r[i]=0;
                deltaRr[i]=0;

                deltaSp[i]=0;
                deltaE1p[i]=0;
                deltaE2p[i]=0;
                deltaI1p[i]=0;
                deltaI2p[i]=0;
                deltaRp[i]=0;

                /*rate of depletion of susceptible*/
                for(j=0;j<NAG;j++)
                    deltaS[i]+=transmission_regular[i*7+j]*(I1[j]+I2[j]+I1r[j]+I2r[j]+I1p[j]+I2p[j]);

                deltaSr[i]=deltaS[i];
                deltaSp[i]=deltaS[i];
                deltaS[i]*=-S[i];
                deltaSr[i]*=-Sr[i];
                deltaSp[i]*=-Sp[i];

                /*rate of passing between states of infection*/
                deltaE1[i]=-deltaS[i]-a1*E1[i];
                deltaE1r[i]=-deltaSr[i]-a1*E1r[i];
                deltaE1p[i]=-deltaSp[i]-a1*E1p[i];
                deltaE2[i]=a1*E1[i]-a2*E2[i];
                deltaE2r[i]=a1*E1r[i]-a2*E2r[i];
                deltaE2p[i]=a1*E1p[i]-a2*E2p[i];
                deltaI1[i]=a2*E2[i]-g1*I1[i];
                deltaI1r[i]=a2*E2r[i]-g1*I1r[i];
                deltaI1p[i]=a2*E2p[i]-g1*I1p[i];
                deltaI2[i]=g1*I1[i]-g2*I2[i];
                deltaI2r[i]=g1*I1r[i]-g2*I2r[i];
                deltaI2p[i]=g1*I1p[i]-g2*I2p[i];
                deltaR[i]=g2*I2[i];
                deltaRr[i]=g2*I2r[i];
                deltaRp[i]=g2*I2p[i];
            }

            /*update the different classes*/
            for(i=0;i<NAG;i++)
            {
                S[i]+=h_step*deltaS[i];
                Sr[i]+=h_step*deltaSr[i];
                Sp[i]+=h_step*deltaSp[i];
                E1[i]+=h_step*deltaE1[i];
                E1r[i]+=h_step*deltaE1r[i];
                E1p[i]+=h_step*deltaE1p[i];
                E2[i]+=h_step*deltaE2[i];
                E2r[i]+=h_step*deltaE2r[i];
                E2p[i]+=h_step*deltaE2p[i];
                I1[i]+=h_step*deltaI1[i];
                I1r[i]+=h_step*deltaI1r[i];
                I1p[i]+=h_step*deltaI1p[i];
                I2[i]+=h_step*deltaI2[i];
                I2r[i]+=h_step*deltaI2r[i];
                I2p[i]+=h_step*deltaI2p[i];
                R[i]+=h_step*deltaR[i];
                Rr[i]+=h_step*deltaRr[i];
                Rp[i]+=h_step*deltaRp[i];
                total_of_new_cases_per_day[i]+=a2*E2[i]*h_step;
                total_of_new_cases_per_day_r[i]+=a2*E2r[i]*h_step;
                total_of_new_cases_per_day_p[i]+=a2*E2p[i]*h_step;
            }

            if((((int)(t*step_rate))%step_rate)==step_rate/2)
            {
                for(i=0;i<NAG;i++)
                {
                    result[(int)t*3*NAG+i]=total_of_new_cases_per_day[i];
                    result[(int)t*3*NAG+NAG+i]=total_of_new_cases_per_day_r[i];
                    result[(int)t*3*NAG+2*NAG+i]=total_of_new_cases_per_day_p[i];
                    total_of_new_cases_per_day[i]=0;
                    total_of_new_cases_per_day_r[i]=0;
                    total_of_new_cases_per_day_p[i]=0;
                }

            }
        }
    }

    void days_to_weeks_5AG(double *result_days, double *result_weeks)
    {
        int i,j;

        /*initialisation*/
        for(i=0; i<length_weeks*5; i++)
            result_weeks[i]=0;

        for(i=0; i<length_weeks; i++)
            for(j=0;j<7;j++)
            {
                result_weeks[i*5]+=result_days[(7*i+j)*NAG*3]+result_days[(7*i+j)*NAG*3+1]+result_days[(7*i+j)*NAG*3+7]+result_days[(7*i+j)*NAG*3+8]+result_days[(7*i+j)*NAG*3+14]+result_days[(7*i+j)*NAG*3+15];
                result_weeks[i*5+1]+=result_days[(7*i+j)*NAG*3+2]+result_days[(7*i+j)*NAG*3+9]+result_days[(7*i+j)*NAG*3+16];
                result_weeks[i*5+2]+=result_days[(7*i+j)*NAG*3+3]+result_days[(7*i+j)*NAG*3+4]+result_days[(7*i+j)*NAG*3+10]+result_days[(7*i+j)*NAG*3+11]+result_days[(7*i+j)*NAG*3+17]+result_days[(7*i+j)*NAG*3+18];
                result_weeks[i*5+3]+=result_days[(7*i+j)*NAG*3+5]+result_days[(7*i+j)*NAG*3+12]+result_days[(7*i+j)*NAG*3+19];
                result_weeks[i*5+4]+=result_days[(7*i+j)*NAG*3+6]+result_days[(7*i+j)*NAG*3+13]+result_days[(7*i+j)*NAG*3+20];
            }
    }

    double log_likelihood_hyper_poisson(const std::vector<double> &eps, double psi, double * result_simu,
            Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
            Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
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
            for(week=0;week<52;week++)
            {
                pop_mon=mon_pop(week,i);
                Z_in_mon=(int)round(result_simu[week*5+i]*pop_mon/pop_5AG_RCGP[i]);
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
