/**
 * @file basic.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 
 * @version 0.1
 * @date 2023-04-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef CLPOLY_PARSE_BASIC_HH
#define CLPOLY_PARSE_BASIC_HH

#include <string>
#include <set>

namespace clpoly{
    namespace parse{
    
 

    inline bool is_blank(char c)
    {
        if (c==' ' || c=='\n' || c=='\r')
        {
            return true;
        } 
        return false;
    }

    inline bool is_parentheses(char c)
    {
        if (c=='(' || c==')')
        {
            return true;
        } 
        return false;
    }
    inline bool is_note(char c)
    {
        return (c==';');
    }
    inline bool is_quoted(char c)
    {
        return (c=='|');
    }
    inline bool is_number(char c)
    {
        return c<='9' && c>='0'; 
    }
    inline bool is_characters(char c)
    {
        static const std::set<char> chars={
            '~', '!', '@', '$', '%', '^', '&', '*', '_', '-', '+', '=', '<', '>', '.' ,'?', '/'
        };
        return chars.find(c)!=chars.end();
    }

    }
}
#endif