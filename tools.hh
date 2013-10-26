#ifndef TOOLS
#define TOOLS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

class Qstring: public std::string{
public:
    Qstring& operator= (std::string &s){
        this->resize(s.size());
        for (size_t i = 0; i< s.size(); i++) (*this)[i] = s[i];
    }
    friend std::istream& operator>> (std::istream &is, Qstring &other){
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

struct Ncconfig{
    std::string fname;
    long current;
};

/* exit when file doesn't exist */
void no_exist(std::string file){
    std::cerr << "Cannot open " << file << " ...\n";
    exit(-1);
};

/* locate a file in control file */
void locate(std::ifstream& infile, std::string str){
    std::string line;
    while(!infile.eof() && line.find(str, 0) == std::string::npos)
        getline(infile, line);
    if (infile.eof()){
        std::cerr << "Cannot find enty: " << str << std::endl;
        exit(-1);
    }
};

#endif
