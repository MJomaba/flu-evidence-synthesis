#include "model.hh"

gsl_rng * r;

void one_year_SEIR_with_vaccination(double * result, double *Npop, double * seeding_infectious, double tlatent, double tinfectious, double *s_profile, double *contact_regular, double q, double *vaccination_calendar, double *vaccine_efficacy)
{
    double transmission_regular[NAG2];
    double S[7], E1[7], E2[7], I1[7], I2[7], R[7];
    double Sr[7], E1r[7], E2r[7], I1r[7], I2r[7], Rr[7];
    double Sp[7], E1p[7], E2p[7], I1p[7], I2p[7], Rp[7];
    double deltaS[7], deltaE1[7], deltaE2[7], deltaI1[7], deltaI2[7], deltaR[7];
    double deltaSr[7], deltaE1r[7], deltaE2r[7], deltaI1r[7], deltaI2r[7], deltaRr[7];
    double deltaSp[7], deltaE1p[7], deltaE2p[7], deltaI1p[7], deltaI2p[7], deltaRp[7];
    double SN[7], E1N[7], E2N[7], I1N[7], I2N[7], RN[7];
    double SrN[7], E1rN[7], E2rN[7], I1rN[7], I2rN[7], RrN[7];
    double SpN[7], E1pN[7], E2pN[7], I1pN[7], I2pN[7], RpN[7];
    double deltaSN[7], deltaE1N[7], deltaE2N[7], deltaI1N[7], deltaI2N[7], deltaRN[7];
    double deltaSrN[7], deltaE1rN[7], deltaE2rN[7], deltaI1rN[7], deltaI2rN[7], deltaRrN[7];
    double deltaSpN[7], deltaE1pN[7], deltaE2pN[7], deltaI1pN[7], deltaI2pN[7], deltaRpN[7];
    double total_of_new_cases_per_day[7], total_of_new_cases_per_day_r[7], total_of_new_cases_per_day_p[7];
    double vacc_prov, vacc_prov_p, vacc_prov_r;
    int i, j, cal_time;
    double a1, a2, g1, g2, t /*, surv[7]={0,0,0,0,0,0,0}*/;
    int step_rate;

    step_rate=(int)(1/h_step);

    a1=2/tlatent;
    a2=a1;
    g1=2/tinfectious;
    g2=g1;

    /*initialisation, transmission matrix*/
    for(i=0;i<NAG2;i++)
    {
        transmission_regular[i]=q*contact_regular[i]*s_profile[i/7];
    }

    /*initialisation, S,E,I,R*/
    for(i=0;i<NAG;i++)
    {
        E1[i]=0;
        E1r[i]=0;
        E1p[i]=0;
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
        S[i]=0;
        Sr[i]=0;
        Sp[i]=0;

        E1N[i]=seeding_infectious[i]*Npop[i]/(Npop[i]+Npop[i+NAG]+Npop[i+2*NAG]);
        E1rN[i]=seeding_infectious[i]*Npop[i+NAG]/(Npop[i]+Npop[i+NAG]+Npop[i+2*NAG]);
        E1pN[i]=seeding_infectious[i]*Npop[i+2*NAG]/(Npop[i]+Npop[i+NAG]+Npop[i+2*NAG]);
        E2N[i]=0;
        E2rN[i]=0;
        E2pN[i]=0;
        I1N[i]=0;
        I1rN[i]=0;
        I1pN[i]=0;
        I2N[i]=0;
        I2rN[i]=0;
        I2pN[i]=0;
        RN[i]=0;
        RrN[i]=0;
        RpN[i]=0;
        SN[i]=Npop[i]-E1N[i];
        SrN[i]=Npop[i+NAG]-E1rN[i];
        SpN[i]=Npop[i+2*NAG]-E1pN[i];

        total_of_new_cases_per_day[i]=0;
        total_of_new_cases_per_day_r[i]=0;
        total_of_new_cases_per_day_p[i]=0;
    }

    for(t=0; t<length; t+=h_step)
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

            deltaSN[i]=0;
            deltaE1N[i]=0;
            deltaE2N[i]=0;
            deltaI1N[i]=0;
            deltaI2N[i]=0;
            deltaRN[i]=0;

            deltaSrN[i]=0;
            deltaE1rN[i]=0;
            deltaE2rN[i]=0;
            deltaI1rN[i]=0;
            deltaI2rN[i]=0;
            deltaRrN[i]=0;

            deltaSpN[i]=0;
            deltaE1pN[i]=0;
            deltaE2pN[i]=0;
            deltaI1pN[i]=0;
            deltaI2pN[i]=0;
            deltaRpN[i]=0;

            /*rate of depletion of susceptible*/
            for(j=0;j<NAG;j++)
                deltaS[i]+=transmission_regular[i*7+j]*(I1[j]+I2[j]+I1r[j]+I2r[j]+I1p[j]+I2p[j]+I1N[j]+I2N[j]+I1rN[j]+I2rN[j]+I1pN[j]+I2pN[j]);

            deltaSr[i]=deltaS[i];
            deltaSp[i]=deltaS[i];
            deltaSN[i]=deltaS[i];
            deltaSrN[i]=deltaS[i];
            deltaSpN[i]=deltaS[i];
            deltaS[i]*=-S[i];
            deltaSr[i]*=-Sr[i];
            deltaSp[i]*=-Sp[i];
            deltaSN[i]*=-SN[i];
            deltaSrN[i]*=-SrN[i];
            deltaSpN[i]*=-SpN[i];

            /*rate of passing between states of infection*/
            deltaE1[i]=-deltaS[i]-a1*E1[i];
            deltaE1r[i]=-deltaSr[i]-a1*E1r[i];
            deltaE1p[i]=-deltaSp[i]-a1*E1p[i];

            deltaE1N[i]=-deltaSN[i]-a1*E1N[i];
            deltaE1rN[i]=-deltaSrN[i]-a1*E1rN[i];
            deltaE1pN[i]=-deltaSpN[i]-a1*E1pN[i];

            deltaE2[i]=a1*E1[i]-a2*E2[i];
            deltaE2r[i]=a1*E1r[i]-a2*E2r[i];
            deltaE2p[i]=a1*E1p[i]-a2*E2p[i];

            deltaE2N[i]=a1*E1N[i]-a2*E2N[i];
            deltaE2rN[i]=a1*E1rN[i]-a2*E2rN[i];
            deltaE2pN[i]=a1*E1pN[i]-a2*E2pN[i];

            deltaI1[i]=a2*E2[i]-g1*I1[i];
            deltaI1r[i]=a2*E2r[i]-g1*I1r[i];
            deltaI1p[i]=a2*E2p[i]-g1*I1p[i];

            deltaI1N[i]=a2*E2N[i]-g1*I1N[i];
            deltaI1rN[i]=a2*E2rN[i]-g1*I1rN[i];
            deltaI1pN[i]=a2*E2pN[i]-g1*I1pN[i];

            deltaI2[i]=g1*I1[i]-g2*I2[i];
            deltaI2r[i]=g1*I1r[i]-g2*I2r[i];
            deltaI2p[i]=g1*I1p[i]-g2*I2p[i];

            deltaI2N[i]=g1*I1N[i]-g2*I2N[i];
            deltaI2rN[i]=g1*I1rN[i]-g2*I2rN[i];
            deltaI2pN[i]=g1*I1pN[i]-g2*I2pN[i];

            deltaR[i]=g2*I2[i];
            deltaRr[i]=g2*I2r[i];
            deltaRp[i]=g2*I2p[i];

            deltaRN[i]=g2*I2N[i];
            deltaRrN[i]=g2*I2rN[i];
            deltaRpN[i]=g2*I2pN[i];
        }

        /*Vaccine bit*/
        if((t>=44)&(t<167))
        for(i=0;i<NAG;i++)
        {
            cal_time=(int)(t)-44;
            vacc_prov=Npop[i]*vaccination_calendar[cal_time*21+i]/(SN[i]+E1N[i]+E2N[i]+I1N[i]+I2N[i]+RN[i]);
            /*surv[i]+=vaccination_calendar[cal_time*21+i];*/
            vacc_prov_r=Npop[i+NAG]*vaccination_calendar[cal_time*21+i+NAG]/(SrN[i]+E1rN[i]+E2rN[i]+I1rN[i]+I2rN[i]+RrN[i]);
            vacc_prov_p=0; /*Npop[i+2*NAG]*vaccination_calendar[cal_time*21+i+2*NAG]/(SpN[i]+E1pN[i]+E2pN[i]+I1pN[i]+I2pN[i]+RpN[i]);*/

            deltaS[i]+=SN[i]*vacc_prov*(1-vaccine_efficacy[i]);
            deltaSr[i]+=SrN[i]*vacc_prov_r*(1-vaccine_efficacy[i]);
            deltaSp[i]+=SpN[i]*vacc_prov_p*(1-vaccine_efficacy[i]);
            deltaSN[i]-=SN[i]*vacc_prov;
            deltaSrN[i]-=SrN[i]*vacc_prov_r;
            deltaSpN[i]-=SpN[i]*vacc_prov_p;

            deltaE1[i]+=E1N[i]*vacc_prov;
            deltaE1r[i]+=E1rN[i]*vacc_prov_r;
            deltaE1p[i]+=E1pN[i]*vacc_prov_p;
            deltaE1N[i]-=E1N[i]*vacc_prov;
            deltaE1rN[i]-=E1rN[i]*vacc_prov_r;
            deltaE1pN[i]-=E1pN[i]*vacc_prov_p;

            deltaE2[i]+=E2N[i]*vacc_prov;
            deltaE2r[i]+=E2rN[i]*vacc_prov_r;
            deltaE2p[i]+=E2pN[i]*vacc_prov_p;
            deltaE2N[i]-=E2N[i]*vacc_prov;
            deltaE2rN[i]-=E2rN[i]*vacc_prov_r;
            deltaE2pN[i]-=E2pN[i]*vacc_prov_p;

            deltaI1[i]+=I1N[i]*vacc_prov;
            deltaI1r[i]+=I1rN[i]*vacc_prov_r;
            deltaI1p[i]+=I1pN[i]*vacc_prov_p;
            deltaI1N[i]-=I1N[i]*vacc_prov;
            deltaI1rN[i]-=I1rN[i]*vacc_prov_r;
            deltaI1pN[i]-=I1pN[i]*vacc_prov_p;

            deltaI2[i]+=I2N[i]*vacc_prov;
            deltaI2r[i]+=I2rN[i]*vacc_prov_r;
            deltaI2p[i]+=I2pN[i]*vacc_prov_p;
            deltaI2N[i]-=I2N[i]*vacc_prov;
            deltaI2rN[i]-=I2rN[i]*vacc_prov_r;
            deltaI2pN[i]-=I2pN[i]*vacc_prov_p;

            deltaR[i]+=RN[i]*vacc_prov+SN[i]*vacc_prov*vaccine_efficacy[i];
            deltaRr[i]+=RrN[i]*vacc_prov_r+SrN[i]*vacc_prov_r*vaccine_efficacy[i];
            deltaRp[i]+=RpN[i]*vacc_prov_p+SpN[i]*vacc_prov_p*vaccine_efficacy[i];
            deltaRN[i]-=RN[i]*vacc_prov;
            deltaRrN[i]-=RrN[i]*vacc_prov_r;
            deltaRpN[i]-=RpN[i]*vacc_prov_p;
        }

        /*update the different classes*/
        for(i=0;i<NAG;i++)
        {
            S[i]+=h_step*deltaS[i];
            Sr[i]+=h_step*deltaSr[i];
            Sp[i]+=h_step*deltaSp[i];

            SN[i]+=h_step*deltaSN[i];
            SrN[i]+=h_step*deltaSrN[i];
            SpN[i]+=h_step*deltaSpN[i];

            E1[i]+=h_step*deltaE1[i];
            E1r[i]+=h_step*deltaE1r[i];
            E1p[i]+=h_step*deltaE1p[i];

            E1N[i]+=h_step*deltaE1N[i];
            E1rN[i]+=h_step*deltaE1rN[i];
            E1pN[i]+=h_step*deltaE1pN[i];

            E2[i]+=h_step*deltaE2[i];
            E2r[i]+=h_step*deltaE2r[i];
            E2p[i]+=h_step*deltaE2p[i];

            E2N[i]+=h_step*deltaE2N[i];
            E2rN[i]+=h_step*deltaE2rN[i];
            E2pN[i]+=h_step*deltaE2pN[i];

            I1[i]+=h_step*deltaI1[i];
            I1r[i]+=h_step*deltaI1r[i];
            I1p[i]+=h_step*deltaI1p[i];

            I1N[i]+=h_step*deltaI1N[i];
            I1rN[i]+=h_step*deltaI1rN[i];
            I1pN[i]+=h_step*deltaI1pN[i];

            I2[i]+=h_step*deltaI2[i];
            I2r[i]+=h_step*deltaI2r[i];
            I2p[i]+=h_step*deltaI2p[i];

            I2N[i]+=h_step*deltaI2N[i];
            I2rN[i]+=h_step*deltaI2rN[i];
            I2pN[i]+=h_step*deltaI2pN[i];

            R[i]+=h_step*deltaR[i];
            Rr[i]+=h_step*deltaRr[i];
            Rp[i]+=h_step*deltaRp[i];

            RN[i]+=h_step*deltaRN[i];
            RrN[i]+=h_step*deltaRrN[i];
            RpN[i]+=h_step*deltaRpN[i];

            total_of_new_cases_per_day[i]+=a2*(E2[i]+E2N[i])*h_step;
            total_of_new_cases_per_day_r[i]+=a2*(E2r[i]+E2rN[i])*h_step;
            total_of_new_cases_per_day_p[i]+=a2*(E2p[i]+E2pN[i])*h_step;
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

void one_year_SEIR_without_vaccination(double * result, double *Npop, double * seeding_infectious, double tlatent, double tinfectious, double *s_profile, double *contact_regular, double q)
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
    for(i=0;i<NAG2;i++)
    {
        transmission_regular[i]=q*contact_regular[i]*s_profile[i/7];
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

    for(t=0; t<length; t+=h_step)
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

double log_likelihood_hyper(double* eps, double * result_simu, int * n_ILI, int * mon_popu, int * n_posi, int * n_sampled, double * pop_5AG_RCGP)
{
    int i, k, week;
    int Z_in_mon, n, m, n_plus, max_m_plus;
    double result, epsilon;
    long double u, likelihood_AG_week;

    result=0.0;
    for(i=0;i<5;i++)
    {
        epsilon=eps[i];
        for(week=0;week<52;week++)
        {
            Z_in_mon=(int)round(result_simu[week*5+i]*mon_popu[week*5+i]/pop_5AG_RCGP[i]);
            n=n_sampled[week*5+i];
            n_plus=n_posi[week*5+i];
            m=n_ILI[week*5+i];

            if(Z_in_mon<n_plus)
                return(-10000);

            /*max_n_plus is the min(Z_in_mon, m-n+n_plus)*/
            if(Z_in_mon>m-n+n_plus)
                max_m_plus=m-n+n_plus;
            else
                max_m_plus=Z_in_mon;

            u=pow(epsilon,n_plus)*pow(1-epsilon,Z_in_mon-n_plus);

            if(n_plus<n)
                for(k=n_plus;k<n;k++)
                    u*=m-k;

            if(n_plus>0)
                for(k=0; k<n_plus;k++)
                    u*=Z_in_mon-k;

            likelihood_AG_week=u;

            if(max_m_plus>n_plus)
                for(k=n_plus+1;k<=max_m_plus;k++)
                {
                    u*=(long double)((m-k-n+n_plus+1)*(Z_in_mon-k+1)*epsilon)/((m-k+1)*(k-n_plus)*(1-epsilon));
                    likelihood_AG_week+=u;
                }
            result+=log(likelihood_AG_week);
        }
    }

    return(result);
}

double log_likelihood_hyper_poisson(double* eps, double psi, double * result_simu, int * n_ILI, int * mon_popu, int * n_posi, int * n_sampled, double * pop_5AG_RCGP, int depth)
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
            pop_mon=mon_popu[week*5+i];
            Z_in_mon=(int)round(result_simu[week*5+i]*pop_mon/pop_5AG_RCGP[i]);
            n=n_sampled[week*5+i];
            n_plus=n_posi[week*5+i];
            m=n_ILI[week*5+i];

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

void save_state(const char *name_file, int number, double tl, double ti, double inf_scale, double q_mat, double *s_profile, double *positivity, double *c_prop, double psi, int *bs_part, double *contact_mat, double *epi_1, double lv, double Accept_rate)
{
    FILE *save_file;
    char full_name[80];
    int i;

    for(i=0;name_file[i]!=0;i++)
        full_name[i]=name_file[i];

    full_name[i]=(char)('0'+(number/1000)%10);
    full_name[i+1]=(char)('0'+(number/100)%10);
    full_name[i+2]=(char)('0'+(number/10)%10);
    full_name[i+3]=(char)('0'+number%10);

    full_name[i+4]='.';
    full_name[i+4+1]='s';
    full_name[i+4+2]='t';
    full_name[i+4+3]='m';
    full_name[i+4+4]=0;

    save_file=write_file(full_name);

    /*save the positivity data*/
    fprintf(save_file,"#Latent and infectious period\n");
    fprintf(save_file,"%e %e\n", tl, ti);
    fprintf(save_file,"#Number of initial infectious starting of the epidemics\n");
    fprintf(save_file,"%e\n", inf_scale);
    fprintf(save_file,"#Q factor for the transmission matrix\n");
    fprintf(save_file,"%e\n", q_mat);
    fprintf(save_file,"#Susceptibility profile\n");
    fprintf(save_file,"%e %e %e %e %e %e %e\n", s_profile[0],s_profile[1],s_profile[2],s_profile[3],s_profile[4],s_profile[5],s_profile[6]);
    fprintf(save_file,"#Positivity\n");
    for(i=0;i<52;i++)
        fprintf(save_file,"%.15e %.15e %.15e %.15e %.15e\n",positivity[i*5],positivity[i*5+1],positivity[i*5+2],positivity[i*5+3],positivity[i*5+4]);
    fprintf(save_file,"#Pick up by GP\n");
    fprintf(save_file,"%e %e %e %e %e\n", c_prop[0],c_prop[1],c_prop[2],c_prop[3],c_prop[4]);
    fprintf(save_file,"#Poisson coefficient for outside transmission\n");
    fprintf(save_file,"%e\n",psi);
    fprintf(save_file,"#Bootstrapped contacts\n");
    for(i=0;i<POLY_PART;i++)
        fprintf(save_file,"%d\n",bs_part[i]);
    fprintf(save_file,"#Contact matrix\n");
    for(i=0;i<7;i++)
        fprintf(save_file,"%e %e %e %e %e %e %e\n",contact_mat[i*7],contact_mat[i*7+1],contact_mat[i*7+2],contact_mat[i*7+3],contact_mat[i*7+4],contact_mat[i*7+5],contact_mat[i*7+6]);
    fprintf(save_file,"#Associated epidemic with vaccine\n");
    for(i=0;i<52;i++)
        fprintf(save_file,"%f %f %f %f %f\n",epi_1[5*i],epi_1[5*i+1],epi_1[5*i+2],epi_1[5*i+3],epi_1[5*i+4]);
    fprintf(save_file,"#Current Likelihood\n");
    fprintf(save_file,"%e\n", lv);
    fprintf(save_file,"#Current acceptance rate\n");
    fprintf(save_file,"%e\n", Accept_rate);


    fclose(save_file);
}

void save_scenarii( FILE *Scen1FS, FILE *Scen2FS, int number, double *pop_vec,  double *prop_init_inf,  double tl, double ti, double q_mat, double *s_profile, double *contact_mat, int n_scenarii, double **vaccine_cal,double **vaccine_efficacy_year, char *path, int * First_write)
{
    double result_simu[7644];
    double FinalSize[21];
    int i, j, scen;
    FILE *ScenRest;
    char temp_string[3], filepath[80];

    /*scenario 1*/
    one_year_SEIR_with_vaccination(result_simu, pop_vec, prop_init_inf, tl, ti, s_profile, contact_mat, q_mat, vaccine_cal[0], vaccine_efficacy_year[0]);
	for(j=0; j<21; j++)
	{
		FinalSize[j]=0.0;
	}
	for(i=0;i<364;i++)
    {
        for(j=0; j<21; j++)
        {
            FinalSize[j]+=result_simu[21*i+j];
        }
    }
    for(j=0; j<21; j++)
    {
        fprintf(Scen1FS,"%f ",FinalSize[j]);
    }
    fprintf(Scen1FS,"\n");

    /*scenario 2*/
    one_year_SEIR_without_vaccination(result_simu, pop_vec, prop_init_inf, tl, ti, s_profile, contact_mat, q_mat);
    for(j=0; j<21; j++)
	{
		FinalSize[j]=0.0;
	}
	for(i=0;i<364;i++)
    {
        for(j=0; j<21; j++)
        {
            FinalSize[j]+=result_simu[21*i+j];
        }
    }
    for(j=0; j<21; j++)
    {
        fprintf(Scen2FS,"%f ",FinalSize[j]);
    }
    fprintf(Scen2FS,"\n");

    /*the rest of the scenarios*/
    for(scen=0;scen<n_scenarii;scen++)
    {
        strcpy(filepath,path);
        sprintf(temp_string,"%d",scen);
        strcat(filepath,"scenarii/Scenario_");
        strcat(filepath,temp_string);
        strcat(filepath,"_final_size.txt");
        if(*First_write==1)
        {
            ScenRest=write_file(filepath);
        }
        else
        {
            ScenRest=append_file(filepath);
        }


        one_year_SEIR_with_vaccination(result_simu, pop_vec, prop_init_inf, tl, ti, s_profile, contact_mat, q_mat, vaccine_cal[scen+1], vaccine_efficacy_year[scen+1]);
        for(j=0; j<21; j++)
        {
            FinalSize[j]=0.0;
        }
        for(i=0;i<364;i++)
        {
            for(j=0; j<21; j++)
            {
                FinalSize[j]+=result_simu[21*i+j];
            }
        }
        for(j=0; j<21; j++)
        {
            fprintf(ScenRest,"%f ",FinalSize[j]);
        }
        fprintf(ScenRest,"\n");
        fclose(ScenRest);
    }
    *First_write=0;
}

void proposal_haario(parameter_set * current, parameter_set * proposed, double * chol_de, double * chol_ini, int n, double beta)
{
    double normal_draw[9], normal_add_draw[9], correlated_draw[9], correlated_fix[9];
    double unif1, unif2;
    int i, j, valid_flag;
    double un_moins_beta;

    un_moins_beta=1-beta;

    valid_flag=0;
    do
    {
        /*drawing of the 9 N(0,1) samples using Box-Muller*/
        for(i=0;i<4;i++)
        {
            unif1=gsl_rng_uniform (r);
            unif2=gsl_rng_uniform (r);
            normal_draw[i*2]=2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2); /*3 = sqrt(9)*/
            normal_draw[i*2+1]=2.38/3*sqrt(-2*log(unif1))*cos(twopi*unif2);
        }

        unif1=gsl_rng_uniform (r);
        unif2=gsl_rng_uniform (r);
        normal_draw[8]=2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2);
        normal_add_draw[8]=sqrt(-2*log(unif1))*cos(twopi*unif2);

        /*drawing of the 9 N(0,1) samples using Box-Muller*/
        for(i=0;i<4;i++)
        {
            unif1=gsl_rng_uniform (r);
            unif2=gsl_rng_uniform (r);
            normal_add_draw[i*2]=sqrt(-2*log(unif1))*sin(twopi*unif2);
            normal_add_draw[i*2+1]=sqrt(-2*log(unif1))*cos(twopi*unif2);
        }

        /*transforming the numbers generated with the Cholesky mat to get the correlated samples*/
        for(i=0;i<9;i++)
            correlated_draw[i]=0;

        for(i=0;i<9;i++)
            for(j=0;j<=i;j++)
                correlated_draw[i]+=chol_de[i*9+j]*normal_draw[j];

        for(i=0;i<9;i++)
            correlated_fix[i]=0;

        for(i=0;i<9;i++)
            for(j=0;j<=i;j++)
                correlated_fix[i]+=chol_ini[i*9+j]*normal_add_draw[j];

        /*new proposed values*/
        proposed->epsilon[0]=current->epsilon[0]+un_moins_beta*correlated_draw[0]+beta*correlated_fix[0];
        proposed->epsilon[1]=proposed->epsilon[0];
        proposed->epsilon[2]=current->epsilon[2]+un_moins_beta*correlated_draw[1]+beta*correlated_fix[1];
        proposed->epsilon[3]=proposed->epsilon[2];
        proposed->epsilon[4]=current->epsilon[4]+un_moins_beta*correlated_draw[2]+beta*correlated_fix[2];

        proposed->psi=current->psi+un_moins_beta*correlated_draw[3]+beta*correlated_fix[3];

        proposed->transmissibility=current->transmissibility+un_moins_beta*correlated_draw[4]+beta*correlated_fix[4];

        proposed->susceptibility[0]=current->susceptibility[0]+un_moins_beta*correlated_draw[5]+beta*correlated_fix[5];
        proposed->susceptibility[1]=proposed->susceptibility[0];
        proposed->susceptibility[2]=proposed->susceptibility[0];
        proposed->susceptibility[3]=current->susceptibility[3]+un_moins_beta*correlated_draw[6]+beta*correlated_fix[6];
        proposed->susceptibility[4]=proposed->susceptibility[3];
        proposed->susceptibility[5]=proposed->susceptibility[3];
        proposed->susceptibility[6]=current->susceptibility[6]+un_moins_beta*correlated_draw[7]+beta*correlated_fix[7];

        proposed->init_pop=current->init_pop+un_moins_beta*correlated_draw[8]+beta*normal_add_draw[8];

        /*checking that the generating values are ok i.e. between 0 and 1 if probabilities*/
        if(proposed->epsilon[0] > 0)
            if(proposed->epsilon[0] < 1)
                if(proposed->epsilon[2] > 0)
                    if(proposed->epsilon[2] < 1)
                        if(proposed->epsilon[4] > 0)
                            if(proposed->epsilon[4] < 1)
                                if(proposed->psi >= 0)
                                    if(proposed->psi <= 1)
                                        if(proposed->transmissibility >= 0)
                                            if(proposed->transmissibility <= 1)
                                                if(proposed->susceptibility[0] >= 0)
                                                    if(proposed->susceptibility[0] <= 1)
                                                        if(proposed->susceptibility[3] >= 0)
                                                            if(proposed->susceptibility[3] <= 1)
                                                                if(proposed->susceptibility[6] >= 0)
                                                                    if(proposed->susceptibility[6] <= 1)
                                                                        if(proposed->init_pop<5)
                                                                            valid_flag=1;
    }
    while(valid_flag==0);

}

void proposal_haario_adapt_scale(parameter_set * current, parameter_set * proposed, double * chol_de, double * chol_ini, int n, double beta, double adapt_scale)
{
    double normal_draw[9], normal_add_draw[9], correlated_draw[9], correlated_fix[9];
    double unif1, unif2;
    int i, j, valid_flag;
    double un_moins_beta;

    un_moins_beta=1-beta;

    valid_flag=0;
    do
    {
        /*drawing of the 9 N(0,1) samples using Box-Muller*/
        for(i=0;i<4;i++)
        {
            unif1=gsl_rng_uniform (r);
            unif2=gsl_rng_uniform (r);
            normal_draw[i*2]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2); /*3 = sqrt(9)*/
            normal_draw[i*2+1]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*cos(twopi*unif2);
        }

        unif1=gsl_rng_uniform (r);
        unif2=gsl_rng_uniform (r);
        normal_draw[8]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2);
        normal_add_draw[8]=sqrt(-2*log(unif1))*cos(twopi*unif2);

        /*drawing of the 9 N(0,1) samples using Box-Muller*/
        for(i=0;i<4;i++)
        {
            unif1=gsl_rng_uniform (r);
            unif2=gsl_rng_uniform (r);
            normal_add_draw[i*2]=sqrt(-2*log(unif1))*sin(twopi*unif2);
            normal_add_draw[i*2+1]=sqrt(-2*log(unif1))*cos(twopi*unif2);
        }

        /*transforming the numbers generated with the Cholesky mat to get the correlated samples*/
        for(i=0;i<9;i++)
            correlated_draw[i]=0;

        for(i=0;i<9;i++)
            for(j=0;j<=i;j++)
                correlated_draw[i]+=chol_de[i*9+j]*normal_draw[j];

        for(i=0;i<9;i++)
            correlated_fix[i]=0;

        for(i=0;i<9;i++)
            for(j=0;j<=i;j++)
                correlated_fix[i]+=chol_ini[i*9+j]*normal_add_draw[j];

        /*new proposed values*/
        proposed->epsilon[0]=current->epsilon[0]+un_moins_beta*correlated_draw[0]+beta*correlated_fix[0];
        proposed->epsilon[1]=proposed->epsilon[0];
        proposed->epsilon[2]=current->epsilon[2]+un_moins_beta*correlated_draw[1]+beta*correlated_fix[1];
        proposed->epsilon[3]=proposed->epsilon[2];
        proposed->epsilon[4]=current->epsilon[4]+un_moins_beta*correlated_draw[2]+beta*correlated_fix[2];

        proposed->psi=current->psi+un_moins_beta*correlated_draw[3]+beta*correlated_fix[3];

        proposed->transmissibility=current->transmissibility+un_moins_beta*correlated_draw[4]+beta*correlated_fix[4];

        proposed->susceptibility[0]=current->susceptibility[0]+un_moins_beta*correlated_draw[5]+beta*correlated_fix[5];
        proposed->susceptibility[1]=proposed->susceptibility[0];
        proposed->susceptibility[2]=proposed->susceptibility[0];
        proposed->susceptibility[3]=current->susceptibility[3]+un_moins_beta*correlated_draw[6]+beta*correlated_fix[6];
        proposed->susceptibility[4]=proposed->susceptibility[3];
        proposed->susceptibility[5]=proposed->susceptibility[3];
        proposed->susceptibility[6]=current->susceptibility[6]+un_moins_beta*correlated_draw[7]+beta*correlated_fix[7];

        proposed->init_pop=current->init_pop+un_moins_beta*correlated_draw[8]+beta*normal_add_draw[8];

        /*checking that the generating values are ok i.e. between 0 and 1 if probabilities*/
        if(proposed->epsilon[0] > 0)
            if(proposed->epsilon[0] < 1)
                if(proposed->epsilon[2] > 0)
                    if(proposed->epsilon[2] < 1)
                        if(proposed->epsilon[4] > 0)
                            if(proposed->epsilon[4] < 1)
                                if(proposed->psi >= 0)
                                    if(proposed->psi <= 1)
                                        if(proposed->transmissibility >= 0)
                                            if(proposed->transmissibility <= 1)
                                                if(proposed->susceptibility[0] >= 0)
                                                    if(proposed->susceptibility[0] <= 1)
                                                        if(proposed->susceptibility[3] >= 0)
                                                            if(proposed->susceptibility[3] <= 1)
                                                                if(proposed->susceptibility[6] >= 0)
                                                                    if(proposed->susceptibility[6] <= 1)
                                                                        if(proposed->init_pop<5)
                                                                            valid_flag=1;
    }
    while(valid_flag==0);

}


void cholevsky(double * A, double * res, int dim)
{
    int i, j, k;
    double sum_L2;

    for(i=0;i<dim;i++)
    {
        for(j=i;j<dim;j++)
        {
            sum_L2=A[i*dim+j];
            for(k=0;k<i;k++)
                sum_L2-=res[i*dim+k]*res[j*dim+k];
            if(i==j)
                res[i*dim+i]=sqrt(sum_L2);
            else
                res[j*dim+i]=sum_L2/res[i*dim+i];
        }
    }

}

void update_sum_corr(double * sum_corr, parameter_set * par)
{
    /*first line*/
    sum_corr[0]+=par->epsilon[0]*par->epsilon[0];
    sum_corr[1]+=par->epsilon[0]*par->epsilon[2];
    sum_corr[2]+=par->epsilon[0]*par->epsilon[4];
    sum_corr[3]+=par->epsilon[0]*par->psi;
    sum_corr[4]+=par->epsilon[0]*par->transmissibility;
    sum_corr[5]+=par->epsilon[0]*par->susceptibility[0];
    sum_corr[6]+=par->epsilon[0]*par->susceptibility[3];
    sum_corr[7]+=par->epsilon[0]*par->susceptibility[6];
    sum_corr[8]+=par->epsilon[0]*par->init_pop;

    /*tranpose of first line*/
    sum_corr[9]=sum_corr[1];
    sum_corr[18]=sum_corr[2];
    sum_corr[27]=sum_corr[3];
    sum_corr[36]=sum_corr[4];
    sum_corr[45]=sum_corr[5];
    sum_corr[54]=sum_corr[6];
    sum_corr[63]=sum_corr[7];
    sum_corr[72]=sum_corr[8];

    /*second line*/
    sum_corr[10]+=par->epsilon[2]*par->epsilon[2];
    sum_corr[11]+=par->epsilon[2]*par->epsilon[4];
    sum_corr[12]+=par->epsilon[2]*par->psi;
    sum_corr[13]+=par->epsilon[2]*par->transmissibility;
    sum_corr[14]+=par->epsilon[2]*par->susceptibility[0];
    sum_corr[15]+=par->epsilon[2]*par->susceptibility[3];
    sum_corr[16]+=par->epsilon[2]*par->susceptibility[6];
    sum_corr[17]+=par->epsilon[2]*par->init_pop;

    /*transpose of second line*/
    sum_corr[19]=sum_corr[11];
    sum_corr[28]=sum_corr[12];
    sum_corr[37]=sum_corr[13];
    sum_corr[46]=sum_corr[14];
    sum_corr[55]=sum_corr[15];
    sum_corr[64]=sum_corr[16];
    sum_corr[73]=sum_corr[17];

    /*third line*/
    sum_corr[20]+=par->epsilon[4]*par->epsilon[4];
    sum_corr[21]+=par->epsilon[4]*par->psi;
    sum_corr[22]+=par->epsilon[4]*par->transmissibility;
    sum_corr[23]+=par->epsilon[4]*par->susceptibility[0];
    sum_corr[24]+=par->epsilon[4]*par->susceptibility[3];
    sum_corr[25]+=par->epsilon[4]*par->susceptibility[6];
    sum_corr[26]+=par->epsilon[4]*par->init_pop;

    /*transpose of third line*/
    sum_corr[29]=sum_corr[21];
    sum_corr[38]=sum_corr[22];
    sum_corr[47]=sum_corr[23];
    sum_corr[56]=sum_corr[24];
    sum_corr[65]=sum_corr[25];
    sum_corr[74]=sum_corr[26];

    /*fourth line*/
    sum_corr[30]+=par->psi*par->psi;
    sum_corr[31]+=par->psi*par->transmissibility;
    sum_corr[32]+=par->psi*par->susceptibility[0];
    sum_corr[33]+=par->psi*par->susceptibility[3];
    sum_corr[34]+=par->psi*par->susceptibility[6];
    sum_corr[35]+=par->psi*par->init_pop;

    /*transpose fourth line*/
    sum_corr[39]=sum_corr[31];
    sum_corr[48]=sum_corr[32];
    sum_corr[57]=sum_corr[33];
    sum_corr[66]=sum_corr[34];
    sum_corr[75]=sum_corr[35];

    /*fifth line*/
    sum_corr[40]+=par->transmissibility*par->transmissibility;
    sum_corr[41]+=par->transmissibility*par->susceptibility[0];
    sum_corr[42]+=par->transmissibility*par->susceptibility[3];
    sum_corr[43]+=par->transmissibility*par->susceptibility[6];
    sum_corr[44]+=par->transmissibility*par->init_pop;

    /*transpose fifth line*/
    sum_corr[49]=sum_corr[41];
    sum_corr[58]=sum_corr[42];
    sum_corr[67]=sum_corr[43];
    sum_corr[76]=sum_corr[44];

    /*sixth line*/
    sum_corr[50]+=par->susceptibility[0]*par->susceptibility[0];
    sum_corr[51]+=par->susceptibility[0]*par->susceptibility[3];
    sum_corr[52]+=par->susceptibility[0]*par->susceptibility[6];
    sum_corr[53]+=par->susceptibility[0]*par->init_pop;

    /*transpose sixth line*/
    sum_corr[59]=sum_corr[51];
    sum_corr[68]=sum_corr[52];
    sum_corr[77]=sum_corr[53];

    /*seventh line*/
    sum_corr[60]+=par->susceptibility[3]*par->susceptibility[3];
    sum_corr[61]+=par->susceptibility[3]*par->susceptibility[6];
    sum_corr[62]+=par->susceptibility[3]*par->init_pop;

    /*transpose seventh line*/
    sum_corr[69]=sum_corr[61];
    sum_corr[78]=sum_corr[62];

    /*eigth line*/
    sum_corr[70]+=par->susceptibility[6]*par->susceptibility[6];
    sum_corr[71]+=par->susceptibility[6]*par->init_pop;

    /*transpose eigth line*/
    sum_corr[79]=sum_corr[71];

    /*ninth line*/
    sum_corr[80]+=par->init_pop*par->init_pop;
}

FILE * read_file( const std::string pathname, const std::string filename )
{
    boost::filesystem::path path = pathname;
    path /= filename;
    if (!boost::filesystem::exists( path )) {
        std::cerr << "File does not exist: " << path << std::endl;
    }
    return fopen(path.c_str(),"r");
}

FILE * write_file( const std::string filename )
{
    boost::filesystem::path filepath = filename;
    boost::filesystem::path path = filename;
    // Create directory if it doesn't exist
    // TODO: There is probably a function to get the path without using remove_filename()
    // If so then the copy is not needed anymore
    boost::filesystem::create_directory( path.remove_filename() );
    FILE * file = fopen( filepath.c_str(), "w+t" );
    return file;
}

FILE * append_file( const std::string filename )
{
    boost::filesystem::path filepath = filename;
    boost::filesystem::path path = filename;
    // Create directory if it doesn't exist
    // TODO: There is probably a function to get the path without using remove_filename()
    // If so then the copy is not needed anymore
    boost::filesystem::create_directory( path.remove_filename() );
    if (!boost::filesystem::exists( filepath )) {
        return fopen(filepath.c_str(), "w+t");
    }
    return fopen(filepath.c_str(), "a");
}
