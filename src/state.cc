#include "state.hh"

State load_state( const std::string &file_path )
{
    auto f_init = read_file( file_path );
    State state;


    /*Reading the init file*/
    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf", &tl, &ti);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &current_par->init_pop);

    /*translate into an initial infected population*/
    for(i=0;i<NAG;i++)
        curr_init_inf[i]=pow(10,current_par->init_pop);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &current_par->transmissibility);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf", &current_par->susceptibility[0],&current_par->susceptibility[1],&current_par->susceptibility[2],&current_par->susceptibility[3],&current_par->susceptibility[4],&current_par->susceptibility[5],&current_par->susceptibility[6]);

    save_fgets(sbuffer, 100, f_init);
    for(i=0;i<52;i++)
    {
        save_fgets(sbuffer, 150, f_init);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf",&p_ij[i*5],&p_ij[i*5+1],&p_ij[i*5+2],&p_ij[i*5+3],&p_ij[i*5+4]);
    }

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf", &current_par->epsilon[0],&current_par->epsilon[1],&current_par->epsilon[2],&current_par->epsilon[3],&current_par->epsilon[4]);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &current_par->psi);

    save_fgets(sbuffer, 100, f_init);

    for(i=0;i<90;i++)
        curr_ni[i]=0;
    curr_nwe=0;
    for(i=0; i<POLY_PART; i++)
    {
        save_fgets(sbuffer, 10, f_init);
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

    return state;
}
