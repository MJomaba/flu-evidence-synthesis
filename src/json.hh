#ifndef FLU_JSON_HH
#define FLU_JSON_HH

#include "client/redef_macros.h"
#include "db/json.h"
#include "bson/bson_stream.hh"

/**
 * \brief Some helper function to make converting JSON easier
 */
namespace json {

    /**
     * \brief Helper function to convert type from json string
     *
     * If the type is a custom (user defined) type it needs to implement:
     *
     * \code{.cpp}
     * friend void operator>>( const mongo::BSONElement &el,
     *        parameter_set &pars );
     * \endcode
     *
     * See flu::state_t for examples
     */
    template<class T> 
        T from_json_string( const std::string &json ) 
        {
            T t;
            mongo::BSONObj json_state = mongo::fromjson( json );
            json_state >> t;
            return t;
        }

    /**
     * \brief Helper function to convert type to json string
     *
     * If the type is a custom (user defined) type it needs to implement:
     *
     * \code{.cpp}
        friend mongo::BSONEmitter &operator<<(
                mongo::BSONEmitter &bbuild, const state_t &state );
     * \endcode
     *
     * See flu::state_t for examples
     */
    template<class T>
        std::string to_json_string( const T &object ) 
        {
            mongo::BSONEmitter bbuild;
            bbuild << object;
            auto bobj = bbuild.obj(); 
            return bobj.jsonString( mongo::Strict, 1 );
        }
};

#endif
