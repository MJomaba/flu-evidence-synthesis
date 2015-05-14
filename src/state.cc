#include "state.hh"

#include "io.hh"

gsl_rng * r;

namespace flu {

state_t load_state( const std::string &file_path, 
        const size_t number_age_groups, const size_t dim_poly_part )
{
    char sbuffer[300];
    auto f_init = read_file( file_path );

    state_t state;

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf", &state.time_latent, &state.time_infectious);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &state.parameters.init_pop);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &state.parameters.transmissibility);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf", &state.parameters.susceptibility[0],&state.parameters.susceptibility[1],&state.parameters.susceptibility[2],&state.parameters.susceptibility[3],&state.parameters.susceptibility[4],&state.parameters.susceptibility[5],&state.parameters.susceptibility[6]);

    save_fgets(sbuffer, 100, f_init);
    for(size_t i=0;i<52;i++)
    {
        save_fgets(sbuffer, 150, f_init);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf",&state.positivity_ij[i*5],&state.positivity_ij[i*5+1],&state.positivity_ij[i*5+2],&state.positivity_ij[i*5+3],&state.positivity_ij[i*5+4]);
    }

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf", &state.parameters.epsilon[0],&state.parameters.epsilon[1],&state.parameters.epsilon[2],&state.parameters.epsilon[3],&state.parameters.epsilon[4]);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
	sscanf(sbuffer,"%lf", &state.parameters.psi);

    save_fgets(sbuffer, 100, f_init);

    for(size_t i=0; i<dim_poly_part; i++)
    {
        int nc;
        save_fgets(sbuffer, 10, f_init);
        sscanf(sbuffer,"%d",&nc);
        state.number_contacts.push_back( nc );
    }

    return state;
}

void save_state(const std::string &file_path, const size_t k, const state_t &state, const int * curr_cnt_number, const double *contact_mat, const double * result_by_week, const double lv, const double Accept_rate) 
{
    FILE *save_file;

    save_file=write_file(file_path +
            (char)('0'+(k/1000)%10) +
            (char)('0'+(k/100)%10) +
            (char)('0'+(k/10)%10) +
            (char)('0'+k%10) +".stm");

    /*save the positivity data*/
    fprintf(save_file,"#Latent and infectious period\n");
    fprintf(save_file,"%e %e\n", state.time_latent, state.time_infectious);
    fprintf(save_file,"#Number of initial infectious starting of the epidemics\n");
    fprintf(save_file,"%e\n", state.parameters.init_pop);
    fprintf(save_file,"#Q factor for the transmission matrix\n");
    fprintf(save_file,"%e\n", state.parameters.transmissibility);
    fprintf(save_file,"#Susceptibility profile\n");
    fprintf(save_file,"%e %e %e %e %e %e %e\n", state.parameters.susceptibility[0],state.parameters.susceptibility[1],state.parameters.susceptibility[2],state.parameters.susceptibility[3],state.parameters.susceptibility[4],state.parameters.susceptibility[5],state.parameters.susceptibility[6]);
    fprintf(save_file,"#Positivity\n");
    for(size_t i=0;i<52;i++)
        fprintf(save_file,"%.15e %.15e %.15e %.15e %.15e\n",state.positivity_ij[i*5],state.positivity_ij[i*5+1],state.positivity_ij[i*5+2],state.positivity_ij[i*5+3],state.positivity_ij[i*5+4]);
    fprintf(save_file,"#Pick up by GP\n");
    fprintf(save_file,"%e %e %e %e %e\n", state.parameters.epsilon[0],state.parameters.epsilon[1],state.parameters.epsilon[2],state.parameters.epsilon[3],state.parameters.epsilon[4]);
    fprintf(save_file,"#Poisson coefficient for outside transmission\n");
    fprintf(save_file,"%e\n",state.parameters.psi);
    fprintf(save_file,"#Bootstrapped contacts\n");
    for(size_t i=0;i<POLY_PART;i++)
        fprintf(save_file,"%d\n",curr_cnt_number[i]);
    fprintf(save_file,"#Contact matrix\n");
    for(size_t i=0;i<7;i++)
        fprintf(save_file,"%e %e %e %e %e %e %e\n",contact_mat[i*7],contact_mat[i*7+1],contact_mat[i*7+2],contact_mat[i*7+3],contact_mat[i*7+4],contact_mat[i*7+5],contact_mat[i*7+6]);
    fprintf(save_file,"#Associated epidemic with vaccine\n");
    for(size_t i=0;i<52;i++)
        fprintf(save_file,"%f %f %f %f %f\n",result_by_week[5*i],result_by_week[5*i+1],result_by_week[5*i+2],result_by_week[5*i+3],result_by_week[5*i+4]);
    fprintf(save_file,"#Current Likelihood\n");
    fprintf(save_file,"%e\n", lv);
    fprintf(save_file,"#Current acceptance rate\n");
    fprintf(save_file,"%e\n", Accept_rate);

    fclose(save_file);
}

};
