#ifndef TOOLS
#define TOOLS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

/* string embraced in quotes */
class qstring: public std::string{
public:
    qstring& operator= (std::string &s){
        this->resize(s.size());
        for (size_t i = 0; i< s.size(); i++) (*this)[i] = s[i];
    }
    friend std::istream& operator>> (std::istream &is, qstring &other){
        char c;
        std::string s;
        while ((char) is.get() != '"' && is.good());
        c = is.get();
        while (c != '"' && is.good()){
            s += c;
            c = is.get();
        }
        other = s;
        return is;
    }
};

/* locate a section in control file */
void sclocate(std::ifstream& infile, std::string str){
    std::string line;
    while(!infile.eof() && line.find(str, 0) == std::string::npos)
        getline(infile, line);
    if (infile.eof()){
        std::cerr << "***** ERROR *****" << std::endl;
        std::cerr << "can not find enty: " << str << std::endl;
        std::cerr << "*****************" << std::endl;
        assert(0);
    }
};

#endif
