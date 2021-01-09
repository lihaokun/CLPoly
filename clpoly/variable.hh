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
#include "basic.hh"
#include <unordered_map>
#include <cassert>
namespace clpoly{

    class variable
    {
        private:
	        static std::vector<std::string> variables;
            static std::unordered_map<std::string,std::size_t> name_map;
            static std::vector<std::size_t> free_serial;
            static std::size_t init_variable(const std::string & variable_name);
            static std::size_t init_variable();
            
        public:
            static variable get_variable(const std::string & variable_name);
            static variable get_variable(std::size_t serial);
            //static variable get_variable(char c) = delete;
            static variable new_variable();
            static variable new_variable(const std::string & variable_name);
            static void del_variable(std::size_t serial);        
            static void del_variable(const std::string & variable_name);        

        protected:
	        std::size_t __serial; 
            explicit variable(std::size_t s):__serial(s)
            {                
                assert(s<variables.size());
                // #ifdef DEBUG
                //     if (serial>=variable::variables.size()){
                //        throw std::invalid_argument("Serial of var not exist.");
                //     }
                // #endif 
            }
        public:
            variable():__serial(0) {}
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

            // constexpr operator bool() const {return __serial;}
            // constexpr operator std::size_t()const {return __serial;}
            // constexpr operator int()const {return __serial;}
            // constexpr operator long()const {return __serial;}
            // constexpr operator uint()const {return __serial;}
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
            inline void swap(variable & v)
            {
                std::swap(this->__serial,v.__serial);
            }
    };
    std::ostream& operator<<  (std::ostream& stream, const std::vector<variable>& v) 
    {
        stream<<"[ ";
        for (auto ptr=v.begin();ptr!=v.end();++ptr)
        {
            if (ptr!=v.begin())
                stream<<" , ";
            stream<<*ptr;

        }
        stream<<" ]";
        return stream;
    }
    std::vector<std::string> variable::variables={""};
    std::unordered_map<std::string,std::size_t> variable::name_map={{"",0}};
    std::vector<std::size_t>  variable::free_serial;
    std::size_t variable::init_variable(const std::string & variable_name)
    {
        if (variable::free_serial.empty())
        {
            variable::name_map[variable_name]=variable::variables.size();
            variable::variables.push_back(variable_name);
            return variable::variables.size()-1;
        }
        else
        {
            std::size_t s=variable::free_serial.back();
            variable::free_serial.pop_back();
            variable::name_map[variable_name]=s;
            variable::variables[s]=variable_name;
            return s;
        }
        
        
    }
    std::size_t variable::init_variable()
    {
        if (variable::free_serial.empty())
        {
            std::string variable_name="#var_"+ std::to_string(variable::variables.size());
            variable::name_map[variable_name]=variable::variables.size();
            variable::variables.push_back(variable_name);
            return variable::variables.size()-1;
        }
        else
        {
            std::size_t s=variable::free_serial.back();
            std::string variable_name="#var_"+ std::to_string(s);
            variable::free_serial.pop_back();
            variable::name_map[variable_name]=s;
            variable::variables[s]=variable_name;
            return s;
        }
        
        
    }
    bool is_user_def_variable_name(const std::string & variable_name)
    {
        return !variable_name.empty() && ( ( 'a'<=variable_name[0]<='z') || ( 'A'<=variable_name[0]<='Z'));
    }
    variable variable::get_variable(const std::string & variable_name)
    {
        auto tmp=variable::name_map.find(variable_name);
        if (tmp==variable::name_map.end())
        {
            assert(is_user_def_variable_name(variable_name));
            return variable(variable::init_variable(variable_name));
        }
        return variable(tmp->second);
    }
    variable variable::new_variable(const std::string & variable_name)
    {
        assert(variable::name_map.find(variable_name)==variable::name_map.end());
        assert(is_user_def_variable_name(variable_name));
        return variable(variable::init_variable(variable_name));
    }

    void variable::del_variable(std::size_t serial)
    {
        if (serial<variable::variables.size() && !variable::variables[serial].empty())
        {
            variable::name_map.erase(variable::variables[serial]);
            variable::variables[serial]="";
            free_serial.push_back(serial);
        }
    }
    void variable::del_variable(const std::string & variable_name)
    {
        auto tmp=variable::name_map.find(variable_name);
        if (tmp!=variable::name_map.end())
        {
            variable::del_variable(tmp->second);
        }
    }


    variable variable::get_variable(std::size_t serial)
    {
       assert(serial<variable::variables.size() && !variable::variables[serial].empty());
       return variable(serial);
    }

    variable variable::new_variable()
    {
        return variable(variable::init_variable());
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
