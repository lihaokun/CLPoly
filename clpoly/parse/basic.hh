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

    std::string get_word()
    {
        // std::stringstream sstr;
        std::string ss;
        char c= 0;
        this->get(c);
        if (this->eof())
            return "";
        while (is_blank(c) || is_note(c))
        {
            while (is_blank(c))
            {
                this->get(c);
                if (this->eof())
                {
                    // std::cout<<""<<std::endl;
                    return "";
                }
            }
            if (is_note(c))
            while (c!='\n')
            {
                this->get(c);
                if (this->eof())
                {
                    // std::cout<<""<<std::endl;
                    return "";
                }
            }
            
        }
        if (is_parentheses(c))
        {
            // std::cout<<c<<std::endl;
            return std::string(1,c);
        }
        if (is_quoted(c))
        {
            while (1)
            {
                c=this->peek();
                if (this->eof())
                    this->get_throw();
                this->get(c);
                if (is_quoted(c))
                    break;
                ss.push_back(c);
            }
            return  ss;
        }
        
        while  (1)
        {
            ss.push_back(c);
            c=this->peek();
            if (this->eof() || is_parentheses(c) || 
                is_note(c) || is_blank(c) || is_quoted(c))
                break;
            this->get(c);
            
        }
        // std::cout<<sstr.str()<<std::endl;
        return  ss;
    }
   
    }
}
#endif