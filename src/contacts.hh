#ifndef FLU_CONTACTS_HH
#define FLU_CONTACTS_HH

#include<array>
#include<vector>
#include<string>

#include<state.hh>

namespace flu {

    // Create a contact, holding its age, and contacts (we, N1, N2, etc)..
    struct contact_t
    {
        size_t id;
        int age, we, N1, N2, N3, N4, N5, N6, N7, AG;
    };

    // Keep ni separate
    struct contacts_t
    {
        std::array<contact_t,POLY_PART> contacts;
        //TODO: What are ni and nwe used for (nwe is not weekend?)
        int ni[90], nwe;
    };

    std::vector<double> load_contact_regular( const std::string& contacts_filepath, const state_t &state, int *age_sizes, int *AG_sizes, const contacts_t &contacts );

};
#endif

