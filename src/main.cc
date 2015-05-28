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
#include "io.hh"
#include "state.hh"
#include "data.hh"
#include "contacts.hh"
#include "vaccine.hh"

using namespace flu;

int main(int argc, char *argv[])
{
    int i, j, k, alea1, alea2, mcmc_chain_length, acceptance, burn_in, thinning;
    contacts::contacts_t prop_c;
    double prop_init_inf[NAG];
    double curr_init_inf[NAG];
    int step_mat, freq_sampling, First_write=1;
    char sbuffer[300];
    double correct_prior, correct_prior_con;
    double alea, p_ac_mat;
    double pop_RCGP[5];
    double result[7644], result_by_week[260]; /*21*52 number of new cases per week*/
    int n_pos[260], n_samples[260], ILI[260], mon_pop[260];
    double lv, prop_likelihood;

    double Accept_rate, past_acceptance;
    double adaptive_scaling, conv_scaling;
    double my_acceptance_rate;
    double init_cov_matrix[81], emp_cov_matrix[81], sum_corr_param_matrix[81], chol_emp_cov[81], chol_ini[81];
    double sum_mean_param[9], sum_check;
    parameter_set proposed_par;
    FILE *log_file;
    FILE *f_pos_sample, *f_n_sample, *f_GP, *f_mon_pop, *Scen1FS, *Scen2FS;
    FILE *f_init_cov, *f_final_cov;
    FILE *f_posterior;

    state_t current_state;

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
    f_pos_sample=read_file(data_path,"positivity.txt");

    /*opens the file with the total number of samples for that strain and season*/
    f_n_sample=read_file(data_path,"n_samples.txt");

    /*opens the file with the RCGP ILI numbers*/
    f_GP=read_file(data_path,"ILI.txt");

    /*opens the file with the size of the monitored population*/
    f_mon_pop=read_file(data_path,"mon_pop.txt");

    /*opens the file with the starting state of the covariance matrix for the proposal*/
    f_init_cov=read_file(data_path,"init_cov_matrix.txt");

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

    /*initialises the matrix which holds the \sum \theta \theta^t*/
    for(i=0; i<dim_par2; i++)
        sum_corr_param_matrix[i]=0;

    /*initialises the vector which holds the \sum \theta*/
    for(i=0; i<dim_par; i++)
        sum_mean_param[i]=0;

    /*load the positivity data*/
    for(i=0;i<52;i++)
        save_fscanf(f_pos_sample,"%d %d %d %d %d",&n_pos[i*5],&n_pos[i*5+1],&n_pos[i*5+2],&n_pos[i*5+3],&n_pos[i*5+4]);

    /*load the number of samples data*/
    for(i=0;i<52;i++)
        save_fscanf(f_n_sample,"%d %d %d %d %d",&n_samples[i*5],&n_samples[i*5+1],&n_samples[i*5+2],&n_samples[i*5+3],&n_samples[i*5+4]);

    /*load the number of GP ILI consultation data*/
    for(i=0;i<52;i++)
        save_fscanf(f_GP,"%d %d %d %d %d",&ILI[i*5],&ILI[i*5+1],&ILI[i*5+2],&ILI[i*5+3],&ILI[i*5+4]);

    /*load the size of the monitored population by week*/
    for(i=0;i<52;i++)
        save_fscanf(f_mon_pop,"%d %d %d %d %d",&mon_pop[i*5],&mon_pop[i*5+1],&mon_pop[i*5+2],&mon_pop[i*5+3],&mon_pop[i*5+4]);


    /*opens the file with the number of positive samples for that strain and season*/
    auto pop_vec = data::load_population( data_path 
            + "age_groups_model.txt" );
    /*pop RCGP*/

    pop_RCGP[0]=pop_vec[0]+pop_vec[1]+pop_vec[7]+pop_vec[8]+pop_vec[14]+pop_vec[15];
    pop_RCGP[1]=pop_vec[2]+pop_vec[9]+pop_vec[16];
    pop_RCGP[2]=pop_vec[3]+pop_vec[4]+pop_vec[10]+pop_vec[11]+pop_vec[17]+pop_vec[18];
    pop_RCGP[3]=pop_vec[5]+pop_vec[12]+pop_vec[19];
    pop_RCGP[4]=pop_vec[6]+pop_vec[13]+pop_vec[20];

    auto c = contacts::load_contacts( data_path + "contacts_for_inference.txt" );

    for(i=0; i<10; i++)
        printf("%d %d %d %d %d %d %d %d %d\n", c.contacts[i].age, c.contacts[i].we, c.contacts[i].N1, c.contacts[i].N2, c.contacts[i].N3, c.contacts[i].N4, c.contacts[i].N5, c.contacts[i].N6, c.contacts[i].N7);
    printf("UK POLYMOD contacts downloaded OK.\n");

    /*definition of the age groups:  0-1 1-4 5-14 15-24 25-44 45-64 65+ */
    auto age_data = data::load_age_data( data_path + "age_sizes.txt" );
    /*printf("Number of w/e days: %d\n",curr_c.nwe);*/

    for(i=0; i<10; i++)
        printf("%lu\n",age_data.age_sizes[i]);
    printf("Age sizes downloaded OK.\n");
    /*end of loading the contacts and sizes of age populations*/

    /*loading parameters regarding vaccination*/
    auto vaccine_programme = vaccine::load_vaccine_programme( 
            data_path+"vaccine_calendar.txt");
    printf("Vaccine calendar loaded.\n");

    /*load the size of the monitored population by week*/
    for(i=0;i<9;i++)
    {
        save_fgets(sbuffer, 300, f_init_cov);
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
    current_state = load_state( data_path + "init_MCMC.txt", 
            NAG, POLY_PART );

    /*translate into an initial infected population*/
    for(i=0;i<NAG;i++)
        curr_init_inf[i]=pow(10,current_state.parameters.init_pop);

    auto curr_c = contacts::shuffle_by_id( c, 
            current_state.number_contacts );

    auto current_contact_regular = 
        contacts::to_symmetric_matrix( curr_c, age_data );

    one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_programme[0] );

    /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
    days_to_weeks_5AG(result,result_by_week);

    /*curr_psi=0.00001;*/
        lv=log_likelihood_hyper_poisson(current_state.parameters.epsilon, current_state.parameters.psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

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

    past_acceptance=acceptance=234;
    freq_sampling=10*thinning;

    conv_scaling=0.001;
    adaptive_scaling=0.3;
    /*mcmc_chain_length=10000;
    acceptance=234;
    freq_sampling=100;
    burn_in=-1;*/

    for(k=1; k<=mcmc_chain_length + burn_in; k++)
    {
        /*update of the variance-covariance matrix and the mean vector*/
        update_sum_corr(sum_corr_param_matrix, &current_state.parameters);
        sum_mean_param[0]+=current_state.parameters.epsilon[0];
        sum_mean_param[1]+=current_state.parameters.epsilon[2];
        sum_mean_param[2]+=current_state.parameters.epsilon[4];
        sum_mean_param[3]+=current_state.parameters.psi;
        sum_mean_param[4]+=current_state.parameters.transmissibility;
        sum_mean_param[5]+=current_state.parameters.susceptibility[0];
        sum_mean_param[6]+=current_state.parameters.susceptibility[3];
        sum_mean_param[7]+=current_state.parameters.susceptibility[6];
        sum_mean_param[8]+=current_state.parameters.init_pop;

        if(k%thinning==0 && k>burn_in)
        {
            f_posterior=append_file(data_path + "posterior.txt");
            fprintf(f_posterior,"%d %e %e %e %e %e %e %e %e %e %e %e\n",k, current_state.parameters.epsilon[0],current_state.parameters.epsilon[2],current_state.parameters.epsilon[4],current_state.parameters.psi,current_state.parameters.transmissibility,current_state.parameters.susceptibility[0],current_state.parameters.susceptibility[3],current_state.parameters.susceptibility[6],current_state.parameters.init_pop,adaptive_scaling,lv);
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
                curr_init_inf[i]=pow(10,current_state.parameters.init_pop);

            one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_programme[0] );
            days_to_weeks_5AG(result,result_by_week);
            /*lv=log_likelihood_hyper_poisson(current_state.parameters.epsilon, current_state.parameters.psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);*/
            Accept_rate=(double)past_acceptance/1000;
    
            for( size_t i = 0; i < POLY_PART; ++i )
                current_state.number_contacts[i] = curr_c.contacts[i].id;
            
            save_state((data_path + "samples/z_hyper").c_str(), k, current_state, current_contact_regular, result_by_week, lv, Accept_rate);
        }

        /*proposal_haario(current_state.parameters,proposed_par,chol_emp_cov,chol_ini,100,0.05);*/
        // TODO: really should make sure proposed_par is equal to current_state.parameters
        // at this point (will result in test failures though)
        proposal_haario_adapt_scale(&current_state.parameters,&proposed_par,chol_emp_cov,chol_ini,100,0.05,adaptive_scaling);

        /*translate into an initial infected population*/
        for(i=0;i<NAG;i++)
            prop_init_inf[i]=pow(10,proposed_par.init_pop);

        /*new proposed contact matrix*/
        /*start from current one*/
        for(i=0;i<POLY_PART;i++)
        {
            prop_c.contacts[i]=curr_c.contacts[i];
        }

        for(i=0;i<90;i++)
            prop_c.ni[i]=curr_c.ni[i];

        prop_c.nwe=curr_c.nwe;

        /*do swap of contacts step_mat times (reduce or increase to change 'distance' of new matrix from current)*/
        if(gsl_rng_uniform (r)<p_ac_mat)
        for(i=0;i<step_mat;i++)
        {
            alea1=gsl_rng_get(r)%POLY_PART;
            alea2=gsl_rng_get(r)%POLY_PART;

            prop_c.ni[prop_c.contacts[alea1].age]--;
            if(prop_c.contacts[alea1].we>0) prop_c.nwe--;

            prop_c.contacts[alea1]=c.contacts[alea2];

            prop_c.ni[prop_c.contacts[alea1].age]++;
            if(prop_c.contacts[alea1].we>0) prop_c.nwe++;
        }
        auto prop_contact_regular = 
            contacts::to_symmetric_matrix( prop_c, age_data );

        /*one_year_SEIR_without_vaccination(result, pop_vec, prop_init_inf, prop_current_state.time_latent, prop_ti, prop_s_profile, prop_contact_regular, prop_q);*/
        one_year_SEIR_with_vaccination(result, pop_vec, prop_init_inf, current_state.time_latent, current_state.time_infectious, proposed_par.susceptibility, prop_contact_regular, proposed_par.transmissibility, vaccine_programme[0] );

        /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
        days_to_weeks_5AG(result,result_by_week);

        /*computes the associated likelihood with the proposed values*/
        prop_likelihood=log_likelihood_hyper_poisson(proposed_par.epsilon, proposed_par.psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

        /*Add the priors*/

        /*correct_prior=0;*/

        // TODO Edwin: reintroduce different priors
        correct_prior=(current_state.parameters.transmissibility-proposed_par.transmissibility)*(current_state.parameters.transmissibility+proposed_par.transmissibility-0.3306366)*650.2099;
        //if(env!=37)
        //{
        //    /*Prior for the transmissibility; year other than 2003/04*/
        //    /*correction for a normal prior with mu=0.1653183 and sd=0.02773053*/
        //    /*prior on q*/
        //    correct_prior=(current_state.parameters.transmissibility-proposed_par.transmissibility)*(current_state.parameters.transmissibility+proposed_par.transmissibility-0.3306366)*650.2099;
        //}

		//if(env==37)
		//{
        //    /*prior on the susceptibility (year 2003/04)*/

        //    /*correction for a normal prior with mu=0.688 and sd=0.083 for the 0-14 */
        //    correct_prior=(current_state.parameters.susceptibility[0]-proposed_par.susceptibility[0])*(current_state.parameters.susceptibility[0]+proposed_par.susceptibility[0]-1.376)*145.1589/2;
        //    /*correction for a normal prior with mu=0.529 and sd=0.122 for the 15-64 */
        //    correct_prior+=(current_state.parameters.susceptibility[3]-proposed_par.susceptibility[3])*(current_state.parameters.susceptibility[3]+proposed_par.susceptibility[3]-1.058)*67.18624/2;
        //    /*correction for a normal prior with mu=0.523 and sd=0.175 for the 65+ */
        //    correct_prior+=(current_state.parameters.susceptibility[6]-proposed_par.susceptibility[6])*(current_state.parameters.susceptibility[6]+proposed_par.susceptibility[6]-1.046)*32.65306/2;
		//}

        /*Prior for the ascertainment probabilities*/

        /*correct for the prior from serology season (lognormal):"0-14" lm=-4.493789, ls=0.2860455*/
        correct_prior_con= log(current_state.parameters.epsilon[0])-log(proposed_par.epsilon[0])+(log(current_state.parameters.epsilon[0])-log(proposed_par.epsilon[0]))*(log(current_state.parameters.epsilon[0])+log(proposed_par.epsilon[0])+8.987578)*6.110824;

        /*correct for the prior from serology season (lognormal):"15-64" lm=-4.117028, ls=0.4751615*/
        correct_prior_con+= log(current_state.parameters.epsilon[2])-log(proposed_par.epsilon[2])+(log(current_state.parameters.epsilon[2])-log(proposed_par.epsilon[2]))*(log(current_state.parameters.epsilon[2])+log(proposed_par.epsilon[2])+8.234056)*2.21456;

        /*correct for the prior from serology season (lognormal):"65+" lm=-2.977965, ls=1.331832*/
        correct_prior_con+= log(current_state.parameters.epsilon[4])-log(proposed_par.epsilon[4])+(log(current_state.parameters.epsilon[4])-log(proposed_par.epsilon[4]))*(log(current_state.parameters.epsilon[4])+log(proposed_par.epsilon[4])+5.95593)*0.2818844;


        /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/

        my_acceptance_rate=exp(prop_likelihood-lv+correct_prior+correct_prior_con);
		if((proposed_par.init_pop<log(0.00001))|(proposed_par.init_pop>log(10)))  my_acceptance_rate=0.0; /* limit the number of initially infected to 10^log(0.00001)> 10^-12 */

        alea=gsl_rng_uniform (r);

        if(alea<my_acceptance_rate) /*with prior*/
        {
            /*update the acceptance rate*/
            acceptance++;
            if(k>=1000)
                adaptive_scaling+=0.766*conv_scaling;

            // TODO: Normally we don't use swap, but just set the proposal to current and then get new proposal. Ask mark about this.
            std::swap( current_state.parameters, proposed_par );

            /*update current likelihood*/
            lv=prop_likelihood;

            /*new proposed contact matrix*/
            /*update*/
            for(i=0;i<POLY_PART;i++)
            {
                curr_c.contacts[i]=prop_c.contacts[i];
            }

            for(i=0;i<90;i++)
                curr_c.ni[i]=prop_c.ni[i];

            curr_c.nwe=prop_c.nwe;

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

    fclose(log_file);
    fclose(f_pos_sample);
    fclose(f_GP);
    fclose(f_mon_pop);
    fclose(f_n_sample);
    fclose(f_init_cov);
    fclose(Scen1FS);
    fclose(Scen2FS);
    /*fclose(f_posterior);*/

    /********************************************************************************************************************************************************
    END Free the allocated memory and close open files
    *********************************************************************************************************************************************************/

    return 0;
}
