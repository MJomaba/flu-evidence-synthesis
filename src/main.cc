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
#include "proposal.hh"

using namespace flu;

int main(int argc, char *argv[])
{
    int i, k, mcmc_chain_length, burn_in, thinning;
    double prop_init_inf[NAG];
    double curr_init_inf[NAG];
    int step_mat, freq_sampling;
    double pop_RCGP[5];
    double result[7644], result_by_week[260]; /*21*52 number of new cases per week*/
    int n_pos[260], n_samples[260], ILI[260], mon_pop[260];
    double lv, prop_likelihood;

    double my_acceptance_rate;
    FILE *log_file;
    FILE *f_pos_sample, *f_n_sample, *f_GP, *f_mon_pop;
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
        ("prior-susceptibility", "Use the prior on susceptibility from 2003/04 data. Otherwise use more general prior on transmisibility" )
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

    /*opens the file to save the posterior*/
    f_posterior=write_file( data_path + "posterior.txt" );
    fprintf(f_posterior,"k epsilon_0_14 epsilon_15_64 epsilon_65_plus psi q sigma_0_14 sigma_15_64 sigma_65_plus init_pop adaptive_scaling likelihood\n");
    fclose(f_posterior);

    /*********************************************************************************************************************************************************************************
    Initialisation bit LOADING DATA
    **********************************************************************************************************************************************************************************/

	r = gsl_rng_alloc (gsl_rng_mt19937); /*marsenne twister */
	gsl_rng_set (r, seed);

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
    auto age_data = data::load_age_data( data_path + "age_sizes.txt" );
    auto pop_vec = data::load_population( age_data.age_group_sizes );
    /*pop RCGP*/

    pop_RCGP[0]=pop_vec[0]+pop_vec[1]+pop_vec[7]+pop_vec[8]+pop_vec[14]+pop_vec[15];
    pop_RCGP[1]=pop_vec[2]+pop_vec[9]+pop_vec[16];
    pop_RCGP[2]=pop_vec[3]+pop_vec[4]+pop_vec[10]+pop_vec[11]+pop_vec[17]+pop_vec[18];
    pop_RCGP[3]=pop_vec[5]+pop_vec[12]+pop_vec[19];
    pop_RCGP[4]=pop_vec[6]+pop_vec[13]+pop_vec[20];

    auto contact_data = contacts::load_contacts( data_path + "contacts_for_inference.txt" );

    for(i=0; i<10; i++)
        fprintf(log_file,"%d %d %d %d %d %d %d %d %d\n", contact_data.contacts[i].age, contact_data.contacts[i].weekend, contact_data.contacts[i].N1, contact_data.contacts[i].N2, contact_data.contacts[i].N3, contact_data.contacts[i].N4, contact_data.contacts[i].N5, contact_data.contacts[i].N6, contact_data.contacts[i].N7);
    fprintf(log_file,"UK POLYMOD contacts downloaded OK.\n");

    /*definition of the age groups:  0-1 1-4 5-14 15-24 25-44 45-64 65+ */
    //auto age_data = data::load_age_data( data_path + "age_sizes.txt" );
    /*printf("Number of w/e days: %d\n",curr_c.nwe);*/

    for(i=0; i<10; i++)
        fprintf(log_file,"%lu\n",age_data.age_sizes[i]);
    fprintf(log_file,"Age sizes downloaded OK.\n");
    /*end of loading the contacts and sizes of age populations*/

    /*loading parameters regarding vaccination*/
    auto vaccine_programme = vaccine::load_vaccine_programme( 
            data_path+"vaccine_calendar.txt");
    fprintf(log_file,"Vaccine calendar loaded.\n");

    /*********************************************************************************************************************************************************************************
    END of Initialisation bit  LOADING DATA
    **********************************************************************************************************************************************************************************/

    /********************************************************************************************************************************************************
    initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    /*Reading the init file*/
    current_state = load_state_json( data_path + "init_MCMC.txt" );

    /*translate into an initial infected population*/
    for(i=0;i<NAG;i++)
        curr_init_inf[i]=pow(10,current_state.parameters.init_pop);


    auto curr_c = contacts::shuffle_by_id( contact_data, 
            current_state.contact_ids );

    auto current_contact_regular = 
        contacts::to_symmetric_matrix( curr_c, age_data );

    one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_programme[0] );

    /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
    days_to_weeks_5AG(result,result_by_week);

    /*curr_psi=0.00001;*/
    lv = log_likelihood_hyper_poisson(current_state.parameters.epsilon, current_state.parameters.psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

    auto proposal_state = proposal::load( data_path+"init_cov_matrix.txt",
            9 );

    fprintf(log_file,"Initial covariance matrix loaded.\n");

    fprintf(log_file,"====%f====\n",lv);

    /********************************************************************************************************************************************************
    END initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    /**************************************************************************************************************************************************************
    Start of the MCMC
    **************************************************************************************************************************************************************/

    step_mat=1;            /*number of contacts exchanged*/
    double p_ac_mat=0.10;          /*prob to redraw matrices*/

    freq_sampling=10*thinning;

    /*mcmc_chain_length=10000;
    acceptance=234;
    freq_sampling=100;
    burn_in=-1;*/
    fclose(log_file); // Write file

    for(k=1; k<=mcmc_chain_length + burn_in; k++)
    {
        /*update of the variance-covariance matrix and the mean vector*/
        proposal_state = proposal::update( std::move( proposal_state ),
                current_state.parameters, k );

        if(k%thinning==0 && k>burn_in)
        {
            f_posterior=append_file(data_path + "posterior.txt");
            fprintf(f_posterior,"%d %e %e %e %e %e %e %e %e %e %e %e\n",k, current_state.parameters.epsilon[0],current_state.parameters.epsilon[2],current_state.parameters.epsilon[4],current_state.parameters.psi,current_state.parameters.transmissibility,current_state.parameters.susceptibility[0],current_state.parameters.susceptibility[3],current_state.parameters.susceptibility[6],current_state.parameters.init_pop,proposal_state.adaptive_scaling,lv);
            fclose(f_posterior);
        }

        if (k%freq_sampling == 0)
        {
            log_file = append_file( data_path + "final.log" );
            fprintf(log_file,"[%d]",(k-burn_in)/freq_sampling);
            fclose(log_file); // Write file
        }

        /*generates a sample of current state and writes to disk*/
        if((k%freq_sampling==0)&&(k>burn_in))
        {
            /*calculate the current initial infected population*/
            for(i=0;i<NAG;i++)
                curr_init_inf[i]=pow(10,current_state.parameters.init_pop);

            one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_programme[0] );
            days_to_weeks_5AG(result,result_by_week);
            /*lv=log_likelihood_hyper_poisson(current_state.parameters.epsilon, current_state.parameters.psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);*/
    
            for( size_t i = 0; i < POLY_PART; ++i )
                current_state.contact_ids[i] = curr_c.contacts[i].id;
            
            std::string kstring = boost::lexical_cast<std::string>( k );
            while (kstring.size()<4)
                kstring = "0" + kstring;
            save_state_json( current_state, data_path + "samples/z_hyper"
                    + kstring + ".stm" );
            //save_state((data_path + "samples/z_hyper").c_str(), k, current_state, current_contact_regular, result_by_week, lv, Accept_rate);
        }

        /*proposal_haario(current_state.parameters,proposed_par,chol_emp_cov,chol_ini,100,0.05);*/
        auto proposed_par = proposal::haario_adapt_scale(
                current_state.parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,100,0.05, 
                proposal_state.adaptive_scaling, r );

        /*translate into an initial infected population*/
        for(i=0;i<NAG;i++)
            prop_init_inf[i]=pow(10,proposed_par.init_pop);

        auto prop_c = curr_c;

        /*do swap of contacts step_mat times (reduce or increase to change 'distance' of new matrix from current)*/
        if(gsl_rng_uniform (r)<p_ac_mat)
            prop_c = contacts::bootstrap_contacts( std::move(prop_c),
                    contact_data, step_mat, r );

        auto prop_contact_regular = 
            contacts::to_symmetric_matrix( prop_c, age_data );

        /*one_year_SEIR_without_vaccination(result, pop_vec, prop_init_inf, prop_current_state.time_latent, prop_ti, prop_s_profile, prop_contact_regular, prop_q);*/
        one_year_SEIR_with_vaccination(result, pop_vec, prop_init_inf, current_state.time_latent, current_state.time_infectious, proposed_par.susceptibility, prop_contact_regular, proposed_par.transmissibility, vaccine_programme[0] );

        /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
        days_to_weeks_5AG(result,result_by_week);

        /*computes the associated likelihood with the proposed values*/
        prop_likelihood=log_likelihood_hyper_poisson(proposed_par.epsilon, proposed_par.psi, result_by_week, ILI, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

        /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/
        // Make sure accept works with -inf prior
        my_acceptance_rate=exp(prop_likelihood-lv+
                log_prior(proposed_par, current_state.parameters, vm.count("prior-susceptibility") ));

        double alea=gsl_rng_uniform (r);

        if(alea<my_acceptance_rate) /*with prior*/
        {
            /*update the acceptance rate*/
            proposal_state.acceptance++;
            if(k>=1000)
                proposal_state.adaptive_scaling
                    += 0.766*proposal_state.conv_scaling;

            current_state.parameters = proposed_par;

            /*update current likelihood*/
            lv=prop_likelihood;

            /*new proposed contact matrix*/
            /*update*/
            curr_c = prop_c;
            current_contact_regular=prop_contact_regular;
        }
        else /*if reject*/
        {
            if(k>=1000)
                proposal_state.adaptive_scaling
                    -=0.234*proposal_state.conv_scaling;
        }

    }
    /**************************************************************************************************************************************************************
    END of the MCMC
    **************************************************************************************************************************************************************/

    /********************************************************************************************************************************************************
    Free the allocated memory and close open files
    *********************************************************************************************************************************************************/

    fclose(f_pos_sample);
    fclose(f_GP);
    fclose(f_mon_pop);
    fclose(f_n_sample);
    /*fclose(f_posterior);*/

    /********************************************************************************************************************************************************
    END Free the allocated memory and close open files
    *********************************************************************************************************************************************************/

    return 0;
}
