#include "utils.h"
#include <regex.h>

namespace CS {

bool string_match(const std::string &str, const std::string &pattern)
{
    int    status;
    regex_t    re;

    if( regcomp(&re, pattern.c_str(), REG_EXTENDED|REG_NOSUB) != 0 )
        throw "Regular expression syntax error.";

    status = regexec(&re, str.c_str(), (size_t) 0, NULL, 0);
    regfree(&re);

    if (status != 0)
        return false;

    return true;
}

} // namespace CS
