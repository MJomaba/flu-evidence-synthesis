#include "vaccine.hh"

namespace flu {
    namespace vaccine {
        vaccine_t load_vaccine( FILE *file )
        {
            vaccine_t vacc;
            char sbuffer[300];

            sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf",&vacc.efficacy_year[0],&vacc.efficacy_year[1],&vacc.efficacy_year[2],&vacc.efficacy_year[3],&vacc.efficacy_year[4],&vacc.efficacy_year[5],&vacc.efficacy_year[6]);

            save_fgets(sbuffer, 50, file);
            for(size_t j=0;j<123;j++)
            {
                save_fgets(sbuffer, 300, file);
                sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &vacc.calendar[j*21],&vacc.calendar[j*21+1],&vacc.calendar[j*21+2],&vacc.calendar[j*21+3],&vacc.calendar[j*21+4],&vacc.calendar[j*21+5],&vacc.calendar[j*21+6],&vacc.calendar[j*21+7],&vacc.calendar[j*21+8],&vacc.calendar[j*21+9],&vacc.calendar[j*21+10],&vacc.calendar[j*21+11],&vacc.calendar[j*21+12],&vacc.calendar[j*21+13],&vacc.calendar[j*21+14],&vacc.calendar[j*21+15],&vacc.calendar[j*21+16],&vacc.calendar[j*21+17],&vacc.calendar[j*21+18],&vacc.calendar[j*21+19],&vacc.calendar[j*21+20]);
            }
            return vacc;
        }

        std::vector<vaccine_t> load_vaccine_programme( 
                const std::string &path ) 
        {
            auto vacc_programme = read_file( path );
            std::vector<vaccine_t> programme;

            char sbuffer[300];
            save_fgets(sbuffer, 100, vacc_programme);
            save_fgets(sbuffer, 100, vacc_programme);

            int n_scenarii;
            sscanf(sbuffer,"%d",&n_scenarii);

            save_fgets(sbuffer, 100, vacc_programme);
            save_fgets(sbuffer, 100, vacc_programme);
            save_fgets(sbuffer, 100, vacc_programme);
            programme.push_back( load_vaccine( vacc_programme ) );

            for(int i=0;i<n_scenarii;i++)
            {
                save_fgets(sbuffer, 100, vacc_programme);
                save_fgets(sbuffer, 100, vacc_programme);
                save_fgets(sbuffer, 100, vacc_programme);

                programme.push_back( load_vaccine( vacc_programme ) );
            }


            return programme;
        }

    };
};
