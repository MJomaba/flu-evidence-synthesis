#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

#include "model.hh"

int main(int argc, char *argv[])
{
    int i, j, k, age_part,  AG_part, alea1, alea2, mcmc_chain_length, acceptance, nc, burn_in, thinning;
    int c_age[POLY_PART], c_we[POLY_PART], c_N1[POLY_PART], c_N2[POLY_PART], c_N3[POLY_PART], c_N4[POLY_PART], c_N5[POLY_PART], c_N6[POLY_PART], c_N7[POLY_PART], c_AG[POLY_PART], c_ni[90], c_nwe;
    int prop_age[POLY_PART], prop_we[POLY_PART], prop_N1[POLY_PART], prop_N2[POLY_PART], prop_N3[POLY_PART], prop_N4[POLY_PART], prop_N5[POLY_PART], prop_N6[POLY_PART], prop_N7[POLY_PART], prop_AG[POLY_PART], prop_ni[90], prop_nwe, prop_cnt_number[POLY_PART];
    int curr_age[POLY_PART], curr_we[POLY_PART], curr_N1[POLY_PART], curr_N2[POLY_PART], curr_N3[POLY_PART], curr_N4[POLY_PART], curr_N5[POLY_PART], curr_N6[POLY_PART], curr_N7[POLY_PART], curr_AG[POLY_PART], curr_ni[90], curr_nwe, curr_cnt_number[POLY_PART];
    double prop_init_inf[NAG];
    double curr_init_inf[NAG], tl, ti;
    int age_sizes[90], AG_sizes[7], aux, step_mat, freq_sampling, First_write=1;
    char sbuffer[300], path[20]="200001/B/";
    double ww[POLY_PART], mij[49], w_norm[7], cij[49], cij_pro;
    double correct_prior, correct_prior_con;
    double alea, p_ac_mat;
    double current_contact_regular[NAG2], prop_contact_regular[NAG2];
    double pop_vec[21], pop_RCGP[5];
    double result[7644], result_by_week[260]; /*21*52 number of new cases per week*/
    int n_pos[260], n_samples[260], n_scenarii, ILI[260], mon_pop[260];
    double p_ij[260];
    double lv, prop_likelihood;
    double vaccine_efficacy_year[7], vaccine_cal[2583], *VE_pro, *VCAL_pro;
    double *tab_cal[100], *tab_VE[100]; /*so far a maximum of 100 scenarios can be changed of course*/
    double Accept_rate, past_acceptance;
    double adaptive_scaling, conv_scaling;
    double my_acceptance_rate;
    double init_cov_matrix[81], emp_cov_matrix[81], sum_corr_param_matrix[81], chol_emp_cov[81], chol_ini[81];
    double sum_mean_param[9], sum_check;
    parameter_set * pointer_par, * current_par, * proposed_par;
    parameter_set ad_par_1, ad_par_2;
    FILE *log_file, *vacc_programme;
    FILE *f_pos_sample, *f_n_sample, *f_GP, *f_mon_pop, *Scen1FS, *Scen2FS;
    FILE *f_init, *f_init_cov, *f_final_cov;
    FILE *contacts_PM, *pop_sizes, *f_pop_model;
    FILE *f_posterior;

    // Command line options
    namespace po = boost::program_options;
   	po::options_description desc( "Usage: flu-evidence-synthesis --data-path [DIR]" );

    std::string data_path = "./";

    // MCMC variables
    mcmc_chain_length=10000000;
    burn_in=1000000;
    thinning=1000;

    desc.add_options()
        ("help,h", "This message.")
        ("data-path,d", po::value<std::string>( &data_path ), "Path to the data")
        ("chain-length", po::value<int>( &mcmc_chain_length ), "Total mcmc chain length (no. of samples")
        ("burn-in", po::value<int>( &burn_in ), "MCMC burn in period")
        ("thinning", po::value<int>( &thinning ), "Thin the sample by only samples every n steps" )
        ;

    po::variables_map vm;
    po::store( 
            po::command_line_parser( argc, argv ).options(desc).run(),
            vm );

    try {
        po::notify( vm );
    } catch (po::required_option e) {
        std::cout << e.what() << std::endl << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    } 

    data_path = (boost::filesystem::canonical(
                boost::filesystem::complete( data_path ) )).native() + "/";

    /*opens the output log file*/
    log_file= write_file(data_path + "final.log");

    /*opens the file with the number of positive samples for that strain and season*/
    f_pop_model = read_file( data_path, "age_groups_model.txt" );

    /*opens the file with the number of positive samples for that strain and season*/
    f_pos_sample=read_file(data_path,"positivity.txt");

    /*opens the file with the total number of samples for that strain and season*/
    f_n_sample=read_file(data_path,"n_samples.txt");

    /*opens the file with the RCGP ILI numbers*/
    f_GP=read_file(data_path,"ILI.txt");

    /*opens the file with the size of the monitored population*/
    f_mon_pop=read_file(data_path,"mon_pop.txt");

    /*opens the file with the vaccine calendar*/
    vacc_programme=read_file(data_path,"vaccine_calendar.txt");

    /*opens the file with the starting state of the MCMC chain*/
    f_init=read_file(data_path,"init_MCMC.txt");

    /*opens the file with the starting state of the covariance matrix for the proposal*/
    f_init_cov=read_file(data_path,"init_cov_matrix.txt");

    contacts_PM=read_file(data_path,"contacts_for_inference.txt");

    /*opens the file with the different age sizes*/
    pop_sizes=read_file(data_path,"age_sizes.txt");

    /*opens the 1st scenarioFile*/
    Scen1FS=write_file(data_path + "scenarii/Scenario_vaccination_final_size.txt");

    /*opens the 2nd scenarioFile*/
    Scen2FS=write_file(data_path + "scenarii/Scenario_no_vaccination_final_size.txt");

    /*opens the file to save the posterior*/
    f_posterior=write_file( data_path + "posterior.txt" );
    fprintf(f_posterior,"k epsilon_0_14 epsilon_15_64 epsilon_65_plus psi q sigma_0_14 sigma_15_64 sigma_65_plus init_pop adaptive_scaling likelihood\n");
    fclose(f_posterior);

    /*********************************************************************************************************************************************************************************
    Initialisation bit LOADING DATA
    **********************************************************************************************************************************************************************************/

	r = gsl_rng_alloc (gsl_rng_mt19937); /*marsenne twister */
	gsl_rng_set (r, seed);

    /*initialises the addresses of the parameters used by current_par and proposed_par*/
    current_par=&ad_par_1;
    proposed_par=&ad_par_2;

    /*initialises the matrix which holds the \sum \theta \theta^t*/
    for(i=0; i<dim_par2; i++)
        sum_corr_param_matrix[i]=0;

    /*initialises the vector which holds the \sum \theta*/
    for(i=0; i<dim_par; i++)
        sum_mean_param[i]=0;

    /*load the positivity data*/
    for(i=0;i<52;i++)
        fscanf(f_pos_sample,"%d %d %d %d %d",&n_pos[i*5],&n_pos[i*5+1],&n_pos[i*5+2],&n_pos[i*5+3],&n_pos[i*5+4]);

    /*load the number of samples data*/
    for(i=0;i<52;i++)
        fscanf(f_n_sample,"%d %d %d %d %d",&n_samples[i*5],&n_samples[i*5+1],&n_samples[i*5+2],&n_samples[i*5+3],&n_samples[i*5+4]);

    /*load the number of GP ILI consultation data*/
    for(i=0;i<52;i++)
        fscanf(f_GP,"%d %d %d %d %d",&ILI[i*5],&ILI[i*5+1],&ILI[i*5+2],&ILI[i*5+3],&ILI[i*5+4]);

    /*load the size of the monitored population by week*/
    for(i=0;i<52;i++)
        fscanf(f_mon_pop,"%d %d %d %d %d",&mon_pop[i*5],&mon_pop[i*5+1],&mon_pop[i*5+2],&mon_pop[i*5+3],&mon_pop[i*5+4]);

    /*load the size of the age groups for the model that year*/
    fscanf(f_pop_model,"%lf %lf %lf %lf %lf %lf %lf",&pop_vec[0],&pop_vec[1],&pop_vec[2],&pop_vec[3],&pop_vec[4],&pop_vec[5],&pop_vec[6]);

    /*high risk*/
    pop_vec[7]=pop_vec[0]*0.021; /*2.1% in the <1 */
    pop_vec[8]=pop_vec[1]*0.055; /*5.5% in the 1-4 */
    pop_vec[9]=pop_vec[2]*0.098; /*9.8% in the 5-14 */
    pop_vec[10]=pop_vec[3]*0.087; /*8.7% in the 15-24 */
    pop_vec[11]=pop_vec[4]*0.092; /*9.2% in the 25-44 */
    pop_vec[12]=pop_vec[5]*0.183; /*18.3% in the 45-64 */
    pop_vec[13]=pop_vec[6]*0.45; /*45% in the 65+ */

    /*pregnant women (pregnant women not considered in this model*/
    pop_vec[14]=0;
    pop_vec[15]=0;
    pop_vec[16]=0;
    pop_vec[17]=0;
    pop_vec[18]=0;
    pop_vec[19]=0;
    pop_vec[20]=0;

    /*low risk*/
    pop_vec[0]-=pop_vec[7];
    pop_vec[1]-=pop_vec[8];
    pop_vec[2]-=pop_vec[9];
    pop_vec[3]-=pop_vec[10]+pop_vec[17];
    pop_vec[4]-=pop_vec[11]+pop_vec[18];
    pop_vec[5]-=pop_vec[12];
    pop_vec[6]-=pop_vec[13];

    /*pop RCGP*/
    pop_RCGP[0]=pop_vec[0]+pop_vec[1]+pop_vec[7]+pop_vec[8]+pop_vec[14]+pop_vec[15];
    pop_RCGP[1]=pop_vec[2]+pop_vec[9]+pop_vec[16];
    pop_RCGP[2]=pop_vec[3]+pop_vec[4]+pop_vec[10]+pop_vec[11]+pop_vec[17]+pop_vec[18];
    pop_RCGP[3]=pop_vec[5]+pop_vec[12]+pop_vec[19];
    pop_RCGP[4]=pop_vec[6]+pop_vec[13]+pop_vec[20];

    for(i=0; i<90; i++)
        c_ni[i]=0;

    c_nwe=0;

    /*Loading of the participants with their number of contacts from Polymod*/
    for(i=0; i<POLY_PART; i++)
    {
        fscanf(contacts_PM,"%d %d %d %d %d %d %d %d %d", &c_age[i], &c_we[i], &c_N1[i], &c_N2[i], &c_N3[i], &c_N4[i], &c_N5[i], &c_N6[i], &c_N7[i]);
        age_part=c_age[i];
        c_ni[age_part]++;
        if(c_we[i]>0) c_nwe++;
        c_AG[i]=0;
        if(age_part>0) c_AG[i]++;
        if(age_part>4) c_AG[i]++;
        if(age_part>14) c_AG[i]++;
        if(age_part>24) c_AG[i]++;
        if(age_part>44) c_AG[i]++;
        if(age_part>64) c_AG[i]++;
    }

    for(i=0; i<10; i++)
        printf("%d %d %d %d %d %d %d %d %d\n", c_age[i], c_we[i], c_N1[i], c_N2[i], c_N3[i], c_N4[i], c_N5[i], c_N6[i], c_N7[i]);
    printf("UK POLYMOD contacts downloaded OK.\n");

    /*definition of the age groups:  0-1 1-4 5-14 15-24 25-44 45-64 65+ */
    for(i=0; i<7; i++)
        AG_sizes[i]=0;
    for(i=0; i<85; i++)
    {
        fscanf(pop_sizes,"%d",&age_sizes[i]);
        if(i==0)
            AG_sizes[0]=age_sizes[0];
        else
            if(i<5)
                AG_sizes[1]+=age_sizes[i];
            else
                if(i<15)
                    AG_sizes[2]+=age_sizes[i];
                else
                    if(i<25)
                        AG_sizes[3]+=age_sizes[i];
                    else
                        if(i<45)
                            AG_sizes[4]+=age_sizes[i];
                        else
                            if(i<65)
                                AG_sizes[5]+=age_sizes[i];
                            else
                                AG_sizes[6]+=age_sizes[i];
    }
    /*put the remaining (85+) in the older AG*/
    fscanf(pop_sizes,"%d",&aux);
    AG_sizes[6]+=aux;

    /*printf("Number of w/e days: %d\n",curr_nwe);*/

    for(i=0; i<10; i++)
        printf("%d\n",age_sizes[i]);
    printf("Age sizes downloaded OK.\n");
    /*end of loading the contacts and sizes of age populations*/

    /*loading parameters regarding vaccination*/
    fgets(sbuffer, 100, vacc_programme);
    fgets(sbuffer, 100, vacc_programme);
    sscanf(sbuffer,"%d",&n_scenarii);

    fgets(sbuffer, 100, vacc_programme);
    fgets(sbuffer, 100, vacc_programme);
    fgets(sbuffer, 100, vacc_programme);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf",&vaccine_efficacy_year[0],&vaccine_efficacy_year[1],&vaccine_efficacy_year[2],&vaccine_efficacy_year[3],&vaccine_efficacy_year[4],&vaccine_efficacy_year[5],&vaccine_efficacy_year[6]);

    fgets(sbuffer, 50, vacc_programme);
    for(j=0;j<123;j++)
    {
        fgets(sbuffer, 300, vacc_programme);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &vaccine_cal[j*21],&vaccine_cal[j*21+1],&vaccine_cal[j*21+2],&vaccine_cal[j*21+3],&vaccine_cal[j*21+4],&vaccine_cal[j*21+5],&vaccine_cal[j*21+6],&vaccine_cal[j*21+7],&vaccine_cal[j*21+8],&vaccine_cal[j*21+9],&vaccine_cal[j*21+10],&vaccine_cal[j*21+11],&vaccine_cal[j*21+12],&vaccine_cal[j*21+13],&vaccine_cal[j*21+14],&vaccine_cal[j*21+15],&vaccine_cal[j*21+16],&vaccine_cal[j*21+17],&vaccine_cal[j*21+18],&vaccine_cal[j*21+19],&vaccine_cal[j*21+20]);
    }

    tab_cal[0]=vaccine_cal;
    tab_VE[0]=vaccine_efficacy_year;

    for(i=0;i<n_scenarii;i++)
    {
        VE_pro=(double *) malloc(7 * sizeof (double));
        fgets(sbuffer, 100, vacc_programme);
        fgets(sbuffer, 100, vacc_programme);
        fgets(sbuffer, 100, vacc_programme);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf",&VE_pro[0],&VE_pro[1],&VE_pro[2],&VE_pro[3],&VE_pro[4],&VE_pro[5],&VE_pro[6]);
        tab_VE[1+i]=VE_pro;

        VCAL_pro=(double *) malloc(2583 * sizeof (double));
        fgets(sbuffer, 100, vacc_programme);
        for(j=0;j<123;j++)
        {
            fgets(sbuffer, 300, vacc_programme);
            sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &VCAL_pro[j*21],&VCAL_pro[j*21+1],&VCAL_pro[j*21+2],&VCAL_pro[j*21+3],&VCAL_pro[j*21+4],&VCAL_pro[j*21+5],&VCAL_pro[j*21+6],&VCAL_pro[j*21+7],&VCAL_pro[j*21+8],&VCAL_pro[j*21+9],&VCAL_pro[j*21+10],&VCAL_pro[j*21+11],&VCAL_pro[j*21+12],&VCAL_pro[j*21+13],&VCAL_pro[j*21+14],&VCAL_pro[j*21+15],&VCAL_pro[j*21+16],&VCAL_pro[j*21+17],&VCAL_pro[j*21+18],&VCAL_pro[j*21+19],&VCAL_pro[j*21+20]);
        }
        tab_cal[1+i]=VCAL_pro;
    }

    printf("Vaccine calendar loaded.\n");

    /*load the size of the monitored population by week*/
    for(i=0;i<9;i++)
    {
        fgets(sbuffer, 300, f_init_cov);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&init_cov_matrix[i*9],&init_cov_matrix[i*9+1],&init_cov_matrix[i*9+2],&init_cov_matrix[i*9+3],&init_cov_matrix[i*9+4],&init_cov_matrix[i*9+5],&init_cov_matrix[i*9+6],&init_cov_matrix[i*9+7],&init_cov_matrix[i*9+8]);
    }
    printf("Initial covariance matrix loaded.\n");

    /*********************************************************************************************************************************************************************************
    END of Initialisation bit  LOADING DATA
    **********************************************************************************************************************************************************************************/

    /********************************************************************************************************************************************************
    initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    /*Reading the init file*/
    fgets(sbuffer, 100, f_init);
    fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf", &tl, &ti);

    fgets(sbuffer, 100, f_init);
    fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &current_par->init_pop);

    /*translate into an initial infected population*/
    for(i=0;i<NAG;i++)
        curr_init_inf[i]=pow(10,current_par->init_pop);

    fgets(sbuffer, 100, f_init);
    fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &current_par->transmissibility);

    fgets(sbuffer, 100, f_init);
    fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf", &current_par->susceptibility[0],&current_par->susceptibility[1],&current_par->susceptibility[2],&current_par->susceptibility[3],&current_par->susceptibility[4],&current_par->susceptibility[5],&current_par->susceptibility[6]);

    fgets(sbuffer, 100, f_init);
    for(i=0;i<52;i++)
    {
        fgets(sbuffer, 150, f_init);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf",&p_ij[i*5],&p_ij[i*5+1],&p_ij[i*5+2],&p_ij[i*5+3],&p_ij[i*5+4]);
    }

    fgets(sbuffer, 100, f_init);
    fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf", &current_par->epsilon[0],&current_par->epsilon[1],&current_par->epsilon[2],&current_par->epsilon[3],&current_par->epsilon[4]);

    fgets(sbuffer, 100, f_init);
    fgets(sbuffer, 100, f_init);
	sscanf(sbuffer,"%lf", &current_par->psi);

    fgets(sbuffer, 100, f_init);

    for(i=0;i<90;i++)
        curr_ni[i]=0;
    curr_nwe=0;
    for(i=0; i<POLY_PART; i++)
    {
        fgets(sbuffer, 10, f_init);
        sscanf(sbuffer,"%d",&nc);

        age_part=c_age[nc];
        curr_ni[age_part]++;
        if(c_we[nc]>0) curr_nwe++;

        curr_age[i]=c_age[nc];
        curr_AG[i]=c_AG[nc];
        curr_we[i]=c_we[nc];
        curr_N1[i]=c_N1[nc];
        curr_N2[i]=c_N2[nc];
        curr_N3[i]=c_N3[nc];
        curr_N4[i]=c_N4[nc];
        curr_N5[i]=c_N5[nc];
        curr_N6[i]=c_N6[nc];
        curr_N7[i]=c_N7[nc];
        curr_cnt_number[i]=nc;
        prop_cnt_number[i]=nc;
    }

    /*update of the weights*/
    for(i=0;i<7;i++)
        w_norm[i]=0;

    for(i=0;i<49;i++)
        mij[i]=0;

    for(i=0; i<POLY_PART; i++)
    {
        age_part=curr_age[i];
        AG_part=curr_AG[i];
        if(curr_we[i]==0)
            ww[i]=(double)age_sizes[age_part]/curr_ni[age_part]*5/(POLY_PART-curr_nwe);
        else
            ww[i]=(double)age_sizes[age_part]/curr_ni[age_part]*2/curr_nwe;

        w_norm[AG_part]+=ww[i];
        mij[7*AG_part]+=curr_N1[i]*ww[i];
        mij[7*AG_part+1]+=curr_N2[i]*ww[i];
        mij[7*AG_part+2]+=curr_N3[i]*ww[i];
        mij[7*AG_part+3]+=curr_N4[i]*ww[i];
        mij[7*AG_part+4]+=curr_N5[i]*ww[i];
        mij[7*AG_part+5]+=curr_N6[i]*ww[i];
        mij[7*AG_part+6]+=curr_N7[i]*ww[i];
    }

    for(i=0; i<49; i++)
    {
        if(w_norm[i/7]>0)
            mij[i]/=w_norm[i/7];
        cij[i]=mij[i]/AG_sizes[i%7];
    }

    for(i=0; i<7; i++)
    {
        current_contact_regular[i*7+i]=cij[i*7+i];
        for(j=0;j<i;j++)
        {
            cij_pro=(cij[i*7+j]+cij[j*7+i])/2;
            current_contact_regular[i*7+j]=cij_pro;
            current_contact_regular[j*7+i]=cij_pro;
        }
    }

    one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, tl, ti, current_par->susceptibility, current_contact_regular, current_par->transmissibility, vaccine_cal, vaccine_efficacy_year);

    /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
    days_to_weeks_5AG(result,result_by_week);

    /*curr_psi=0.00001;*/
        lv=log_likelihood_hyper_poisson(current_par->epsilon, current_par->psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

    printf("====%f====\n",lv);

    /********************************************************************************************************************************************************
    END initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    /**************************************************************************************************************************************************************
    Start of the MCMC
    **************************************************************************************************************************************************************/

    step_mat=1;            /*number of contacts exchanged*/
    p_ac_mat=0.10;          /*prob to redraw matrices*/

    cholevsky(init_cov_matrix , chol_ini, 9);
    for(i=0;i<81;i++)
        chol_emp_cov[i]=chol_ini[i];

    acceptance=234;
    freq_sampling=10000;

    conv_scaling=0.001;
    adaptive_scaling=0.3;
    /*mcmc_chain_length=10000;
    acceptance=234;
    freq_sampling=100;
    burn_in=-1;*/

    for(k=1; k<=mcmc_chain_length + burn_in; k++)
    {
        /*update of the variance-covariance matrix and the mean vector*/
        update_sum_corr(sum_corr_param_matrix, current_par);
        sum_mean_param[0]+=current_par->epsilon[0];
        sum_mean_param[1]+=current_par->epsilon[2];
        sum_mean_param[2]+=current_par->epsilon[4];
        sum_mean_param[3]+=current_par->psi;
        sum_mean_param[4]+=current_par->transmissibility;
        sum_mean_param[5]+=current_par->susceptibility[0];
        sum_mean_param[6]+=current_par->susceptibility[3];
        sum_mean_param[7]+=current_par->susceptibility[6];
        sum_mean_param[8]+=current_par->init_pop;

        if(k%thinning==0)
        {
            f_posterior=append_file(data_path + "posterior.txt");
            fprintf(f_posterior,"%d %e %e %e %e %e %e %e %e %e %e %e\n",k, current_par->epsilon[0],current_par->epsilon[2],current_par->epsilon[4],current_par->psi,current_par->transmissibility,current_par->susceptibility[0],current_par->susceptibility[3],current_par->susceptibility[6],current_par->init_pop,adaptive_scaling,lv);
            fclose(f_posterior);
        }

        /*adjust variance for MCMC parameters*/
        if(k%1000==0)
        {
            /*accept_adjust=1;
            if(((double)acceptance/1000.0)>0.5) accept_adjust=2;
            if(((double)acceptance/1000.0)<0.15) accept_adjust=0.5;*/

            /*printf("[Sw %d %lf]",acceptance,accept_adjust);*/

            /*update of the adaptive algorithm*/
            for(i=0;i<9;i++)
            {
                emp_cov_matrix[i*9+i]=(sum_corr_param_matrix[i*9+i]-(sum_mean_param[i]*sum_mean_param[i])/k)/(k-1);
                for(j=0;j<i;j++)
                {
                    emp_cov_matrix[i*9+j]=(sum_corr_param_matrix[i*9+j]-(sum_mean_param[i]*sum_mean_param[j])/k)/(k-1);
                    emp_cov_matrix[j*9+i]=emp_cov_matrix[i*9+j];
                }
            }
            cholevsky(emp_cov_matrix,chol_emp_cov,9);

            sum_check=0;
            for(i=0; i<81; i++)
                sum_check+=log(fabs(sum_corr_param_matrix[i]));

            past_acceptance=acceptance;
            acceptance=0;

            conv_scaling/=1.005;
         }

        /*generates a sample*/ // This is inference as far as I can tell
        if((k%freq_sampling==0)&&(k>burn_in))
        {
            printf("[%d]",(k-burn_in)/freq_sampling);

            /*calculate the current initial infected population*/
            for(i=0;i<NAG;i++)
                curr_init_inf[i]=pow(10,current_par->init_pop);
            one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, tl, ti, current_par->susceptibility, current_contact_regular, current_par->transmissibility, vaccine_cal, vaccine_efficacy_year);
            days_to_weeks_5AG(result,result_by_week);
            /*lv=log_likelihood_hyper_poisson(current_par->epsilon, current_par->psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);*/
            Accept_rate=(double)past_acceptance/1000;
            save_state((data_path + "samples/z_hyper").c_str(), (k-burn_in)/freq_sampling, tl, ti, current_par->init_pop, current_par->transmissibility, current_par->susceptibility, p_ij, current_par->epsilon, current_par->psi, curr_cnt_number, current_contact_regular, result_by_week, lv, Accept_rate);
            save_scenarii(Scen1FS, Scen2FS, (k-burn_in)/freq_sampling, pop_vec, curr_init_inf, tl, ti, current_par->transmissibility, current_par->susceptibility, current_contact_regular, n_scenarii, tab_cal, tab_VE, path, &First_write);
        }

        /*proposal_haario(current_par,proposed_par,chol_emp_cov,chol_ini,100,0.05);*/
        proposal_haario_adapt_scale(current_par,proposed_par,chol_emp_cov,chol_ini,100,0.05,adaptive_scaling);

        /*translate into an initial infected population*/
        for(i=0;i<NAG;i++)
            prop_init_inf[i]=pow(10,proposed_par->init_pop);

        /*new proposed contact matrix*/
        /*start from current one*/
        for(i=0;i<POLY_PART;i++)
        {
            prop_age[i]=curr_age[i];
            prop_AG[i]=curr_AG[i];
            prop_we[i]=curr_we[i];
            prop_N1[i]=curr_N1[i];
            prop_N2[i]=curr_N2[i];
            prop_N3[i]=curr_N3[i];
            prop_N4[i]=curr_N4[i];
            prop_N5[i]=curr_N5[i];
            prop_N6[i]=curr_N6[i];
            prop_N7[i]=curr_N7[i];
            prop_cnt_number[i]=curr_cnt_number[i];
        }

        for(i=0;i<90;i++)
            prop_ni[i]=curr_ni[i];

        prop_nwe=curr_nwe;

        /*do swap of contacts step_mat times (reduce or increase to change 'distance' of new matrix from current)*/
        if(gsl_rng_uniform (r)<p_ac_mat)
        for(i=0;i<step_mat;i++)
        {
            alea1=gsl_rng_get(r)%POLY_PART;
            alea2=gsl_rng_get(r)%POLY_PART;

            prop_ni[prop_age[alea1]]--;
            if(prop_we[alea1]>0) prop_nwe--;

            prop_age[alea1]=c_age[alea2];
            prop_we[alea1]=c_we[alea2];
            prop_N1[alea1]=c_N1[alea2];
            prop_N2[alea1]=c_N2[alea2];
            prop_N3[alea1]=c_N3[alea2];
            prop_N4[alea1]=c_N4[alea2];
            prop_N5[alea1]=c_N5[alea2];
            prop_N6[alea1]=c_N6[alea2];
            prop_N7[alea1]=c_N7[alea2];
            prop_AG[alea1]=c_AG[alea2];

            prop_cnt_number[alea1]=alea2;

            prop_ni[prop_age[alea1]]++;
            if(prop_we[alea1]>0) prop_nwe++;
        }

        /*update of the weights*/
        for(i=0;i<7;i++)
            w_norm[i]=0;

        for(i=0;i<49;i++)
            mij[i]=0;

        for(i=0; i<POLY_PART; i++)
        {
            age_part=prop_age[i];
            AG_part=prop_AG[i];
            if(prop_we[i]==0)
                ww[i]=(double)age_sizes[age_part]/prop_ni[age_part]*5/(POLY_PART-prop_nwe);
            else
                ww[i]=(double)age_sizes[age_part]/prop_ni[age_part]*2/prop_nwe;

            w_norm[AG_part]+=ww[i];
            mij[7*AG_part]+=prop_N1[i]*ww[i];
            mij[7*AG_part+1]+=prop_N2[i]*ww[i];
            mij[7*AG_part+2]+=prop_N3[i]*ww[i];
            mij[7*AG_part+3]+=prop_N4[i]*ww[i];
            mij[7*AG_part+4]+=prop_N5[i]*ww[i];
            mij[7*AG_part+5]+=prop_N6[i]*ww[i];
            mij[7*AG_part+6]+=prop_N7[i]*ww[i];
        }

        /*Compute the contact matrix*/
        for(i=0; i<49; i++)
        {
            if(w_norm[i/7]>0)
                mij[i]/=w_norm[i/7];
            cij[i]=mij[i]/AG_sizes[i%7];
        }

        for(i=0; i<7; i++)
        {
            prop_contact_regular[i*7+i]=cij[i*7+i];
            for(j=0;j<i;j++)
            {
                cij_pro=(cij[i*7+j]+cij[j*7+i])/2;
                prop_contact_regular[i*7+j]=cij_pro;
                prop_contact_regular[j*7+i]=cij_pro;
            }
        }

        /*one_year_SEIR_without_vaccination(result, pop_vec, prop_init_inf, prop_tl, prop_ti, prop_s_profile, prop_contact_regular, prop_q);*/
        one_year_SEIR_with_vaccination(result, pop_vec, prop_init_inf, tl, ti, proposed_par->susceptibility, prop_contact_regular, proposed_par->transmissibility, vaccine_cal, vaccine_efficacy_year);

        /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
        days_to_weeks_5AG(result,result_by_week);

        /*computes the associated likelihood with the proposed values*/
        prop_likelihood=log_likelihood_hyper_poisson(proposed_par->epsilon, proposed_par->psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

        /*Add the priors*/

        /*correct_prior=0;*/

        // TODO Edwin: reintroduce different priors
        correct_prior=(current_par->transmissibility-proposed_par->transmissibility)*(current_par->transmissibility+proposed_par->transmissibility-0.3306366)*650.2099;
        //if(env!=37)
        //{
        //    /*Prior for the transmissibility; year other than 2003/04*/
        //    /*correction for a normal prior with mu=0.1653183 and sd=0.02773053*/
        //    /*prior on q*/
        //    correct_prior=(current_par->transmissibility-proposed_par->transmissibility)*(current_par->transmissibility+proposed_par->transmissibility-0.3306366)*650.2099;
        //}

		//if(env==37)
		//{
        //    /*prior on the susceptibility (year 2003/04)*/

        //    /*correction for a normal prior with mu=0.688 and sd=0.083 for the 0-14 */
        //    correct_prior=(current_par->susceptibility[0]-proposed_par->susceptibility[0])*(current_par->susceptibility[0]+proposed_par->susceptibility[0]-1.376)*145.1589/2;
        //    /*correction for a normal prior with mu=0.529 and sd=0.122 for the 15-64 */
        //    correct_prior+=(current_par->susceptibility[3]-proposed_par->susceptibility[3])*(current_par->susceptibility[3]+proposed_par->susceptibility[3]-1.058)*67.18624/2;
        //    /*correction for a normal prior with mu=0.523 and sd=0.175 for the 65+ */
        //    correct_prior+=(current_par->susceptibility[6]-proposed_par->susceptibility[6])*(current_par->susceptibility[6]+proposed_par->susceptibility[6]-1.046)*32.65306/2;
		//}

        /*Prior for the ascertainment probabilities*/

        /*correct for the prior from serology season (lognormal):"0-14" lm=-4.493789, ls=0.2860455*/
        correct_prior_con= log(current_par->epsilon[0])-log(proposed_par->epsilon[0])+(log(current_par->epsilon[0])-log(proposed_par->epsilon[0]))*(log(current_par->epsilon[0])+log(proposed_par->epsilon[0])+8.987578)*6.110824;

        /*correct for the prior from serology season (lognormal):"15-64" lm=-4.117028, ls=0.4751615*/
        correct_prior_con+= log(current_par->epsilon[2])-log(proposed_par->epsilon[2])+(log(current_par->epsilon[2])-log(proposed_par->epsilon[2]))*(log(current_par->epsilon[2])+log(proposed_par->epsilon[2])+8.234056)*2.21456;

        /*correct for the prior from serology season (lognormal):"65+" lm=-2.977965, ls=1.331832*/
        correct_prior_con+= log(current_par->epsilon[4])-log(proposed_par->epsilon[4])+(log(current_par->epsilon[4])-log(proposed_par->epsilon[4]))*(log(current_par->epsilon[4])+log(proposed_par->epsilon[4])+5.95593)*0.2818844;


        /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/

        my_acceptance_rate=exp(prop_likelihood-lv+correct_prior+correct_prior_con);
		if((proposed_par->init_pop<log(0.00001))|(proposed_par->init_pop>log(10)))  my_acceptance_rate=0.0; /* limit the number of initially infected to 10^log(0.00001)> 10^-12 */

        alea=gsl_rng_uniform (r);

        if(alea<my_acceptance_rate) /*with prior*/
        {
            /*update the acceptance rate*/
            acceptance++;
            if(k>=1000)
                adaptive_scaling+=0.766*conv_scaling;

            pointer_par=current_par;
            current_par=proposed_par;
            proposed_par=pointer_par;

            /*update current likelihood*/
            lv=prop_likelihood;

            /*new proposed contact matrix*/
            /*update*/
            for(i=0;i<POLY_PART;i++)
            {
                curr_age[i]=prop_age[i];
                curr_AG[i]=prop_AG[i];
                curr_we[i]=prop_we[i];
                curr_N1[i]=prop_N1[i];
                curr_N2[i]=prop_N2[i];
                curr_N3[i]=prop_N3[i];
                curr_N4[i]=prop_N4[i];
                curr_N5[i]=prop_N5[i];
                curr_N6[i]=prop_N6[i];
                curr_N7[i]=prop_N7[i];
                curr_cnt_number[i]=prop_cnt_number[i];
            }

            for(i=0;i<90;i++)
                curr_ni[i]=prop_ni[i];

            curr_nwe=prop_nwe;

            for(i=0;i<49;i++)
                current_contact_regular[i]=prop_contact_regular[i];
        }
        else /*if reject*/
        {
            if(k>=1000)
                adaptive_scaling-=0.234*conv_scaling;
        }

    }
    /**************************************************************************************************************************************************************
    END of the MCMC
    **************************************************************************************************************************************************************/

    f_final_cov=write_file(data_path + "final_cov.txt");
    for(i=0;i<9;i++)
    {
        for(j=0;j<9;j++)
            fprintf(f_final_cov,"%e ",emp_cov_matrix[i*9+j]);
            fprintf(f_final_cov,"\n");
    }
    fclose(f_final_cov);
    /********************************************************************************************************************************************************
    Free the allocated memory and close open files
    *********************************************************************************************************************************************************/

    for(i=0;i<n_scenarii;i++)
    {
        free(tab_cal[1+i]);
        free(tab_VE[1+i]);
    }

    fclose(f_pop_model);
    fclose(contacts_PM);
    fclose(pop_sizes);
    fclose(log_file);
    fclose(vacc_programme);
    fclose(f_pos_sample);
    fclose(f_GP);
    fclose(f_mon_pop);
    fclose(f_n_sample);
    fclose(f_init);
    fclose(f_init_cov);
    fclose(Scen1FS);
    fclose(Scen2FS);
    /*fclose(f_posterior);*/

    /********************************************************************************************************************************************************
    END Free the allocated memory and close open files
    *********************************************************************************************************************************************************/

    return 0;
}
