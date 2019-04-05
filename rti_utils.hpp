// author: Jungwook Shin @ mgh (jshin13@mgh.harvard.edu)

#ifndef RTI_UTILS_H
#define RTI_UTILS_H

#include <string>
namespace rti{

    //1. string utils
    //from cpluscplus.com
    inline std::string trim_right_copy(const std::string& s,const std::string& delimiters = " \f\n\r\t\v\0\\" )
    {
        return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
    }

    inline std::string trim_left_copy(const std::string& s,const std::string& delimiters = " \f\n\r\t\v\0\\" )
    {
    return s.substr( s.find_first_not_of( delimiters ) );
    }

    inline std::string trim_copy(const std::string& s,const std::string& delimiters = " \f\n\r\t\v\0\\" )
    {
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
    }


}

#endif