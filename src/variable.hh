/*
Module Name:
    variable.hh
Abstract:
    定义类：variable
Author:
    haokun li
Notes:
    variable::new_variable
    variable::get_variable
*/
#ifndef CLPOLY_VARIABLE_HH
#define CLPOLY_VARIABLE_HH
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
namespace clpoly{

    class variable
    {
        private:
	        static std::vector<std::string> variables;
            static std::unordered_map<std::string,std::size_t> name_map;
            static std::size_t init_variable(const std::string & variable_name);
        public:
            static variable get_variable(const std::string & variable_name);
            static variable get_variable(std::size_t serial);
            //static variable get_variable(char c) = delete;
            //static variable new_variable();
            static variable new_variable(const std::string & variable_name);
            

        protected:
	        std::size_t __serial; 
            explicit variable(std::size_t s):__serial(s)
            {
                
                #ifdef DEBUG
                    if (serial>=variable::variables.size()){
                       throw std::invalid_argument("Serial of var not exist.");
                    }
                #endif 
            }
        public:
            variable()
            :__serial(0)
            {
                //std::cout<<"_init_  variable()";
                //init_variable("var_"+ std::to_string(this->__serial));
            }
            // variable(char c)
            // {
            //     this->__serial=get_variable(std::string(1,c)).__serial;
            // }
            variable(const std::string & variable_name)
            {
                this->__serial=get_variable(variable_name).__serial;
            }
            variable(const variable & v)
            {
                this->__serial=v.__serial;
            }
            
            constexpr std::size_t serial() const {return this->__serial;}
            inline const std::string & name() const {return variables[this->__serial];} 
            constexpr bool operator==(const variable & v) const {return this->__serial==v.__serial;}
            constexpr bool operator!=(const variable & v) const {return this->__serial!=v.__serial;}
            constexpr bool operator>(const variable & v) const {return this->__serial>v.__serial;}
            constexpr bool operator<(const variable & v) const {return this->__serial<v.__serial;}
            constexpr bool operator<=(const variable & v) const {return this->__serial<=v.__serial;}
            constexpr bool operator>=(const variable & v) const {return this->__serial>=v.__serial;}
            variable & operator=(const variable & v)
            {
                this->__serial=v.__serial;
                return *this;
            }
            inline const std::string & str() const {return this->name();}
            friend inline std::ostream& operator<<  (std::ostream& stream, const variable& v) {
                return stream<<v.name();
            }
    };
    
    std::vector<std::string> variable::variables={""};
    std::unordered_map<std::string,std::size_t> variable::name_map={{"",0}};
    std::size_t variable::init_variable(const std::string & variable_name)
    {
        #ifdef DEBUG
            if (variable::name_map.find(variable_name)!=variable::name_map.end())
                std::cout<<"Warning: Duplicate variable name.";
        #endif
        variable::name_map[variable_name]=variable::variables.size();
        variable::variables.push_back(variable_name);
        return variable::variables.size()-1;
    }
    variable variable::get_variable(const std::string & variable_name)
    {
        auto tmp=variable::name_map.find(variable_name);
        if (tmp==variable::name_map.end())
        {
            return variable(variable::init_variable(variable_name));
        }
        return variable(tmp->second);
    }
    variable variable::get_variable(std::size_t serial)
    {
        if (serial>=variable::variables.size()){
            throw std::invalid_argument("Serial of var not exist.");
        }
        return variable(serial);
    }
    // variable variable::new_variable()
    // {
    //     return variable(variable::init_variable("var_"+ std::to_string(variable::variables.size())));
    // }
    variable variable::new_variable(const std::string & variable_name)
    {
        return variable(variable::init_variable(variable_name));
    }
}
namespace std{
    template<>
    struct hash<clpoly::variable>
    {
        std::size_t operator()(const clpoly::variable& v) const
        {
            return (hash<size_t>()(v.serial()));
        }
    };
}
#endif
