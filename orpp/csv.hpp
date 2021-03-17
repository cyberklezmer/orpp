#ifndef CSV_HPP
#define CSV_HPP

#include <string>
#include <sstream>
#include "orpp.hpp"

namespace orpp
{

template <char sep=','>
class csv : private std::vector<std::vector<std::string>>
{
    static std::string getstr(std::istream& is)
    {
        char c;
        is.get(c);
        bool quotes = false;
        if(c=='"')
        {
            quotes = true;
            is.get(c);
        }

        std::string r;
        for(;;)
        {
            if((c==sep && !quotes) || is.eof())
                break;
            if(c=='"' && quotes)
            {
               is.get(c);
               if(c!='"')
                  break;
            }
            r.push_back(c);
            is.get(c);
        }
        return r;
    }

    static std::vector<std::vector<std::string>> getvecs(const std::string& fn)
    {
        std::vector<vector<std::string>> res;
        std::ifstream file(fn);
        if(!file.is_open())
        {
            std::ostringstream es;
            es << "cannot open csv file: " << fn ;
            throw exception(es);
        }

        std::string line;
        while (getline(file, line))
        {
            unsigned ls = line.size();
            if(line[ls-1]=='\r')
            {
                if(ls==1)
                    break;
                else
                    line.resize(ls-1);
            }


            std::vector<std::string> row;
            std::istringstream is(line);
            for(;;)
            {
                std::string s = getstr(is);
                row.push_back(s);
                if(is.eof())
                    break;
            }
            res.push_back(row);
        }

        return res;
    }

    static std::string c2excel(unsigned c)
    {
        std::ostringstream os;
        unsigned numchars = 'Z'-'A'+1;
        if(c >= numchars)
            os << static_cast<char>( c / numchars - 1 + 'A' );
        os << static_cast<char>( (c % numchars) + 'A' );
        return os.str();
    }

    static std::string c2excel(unsigned r, unsigned c)
    {
        std::ostringstream os;
        os << c2excel(c) << r+1;
        return os.str();
    }

public:
    unsigned r() const { return (*this).size();}
    unsigned c(unsigned r) const
    {
        if(r>=(*this).size())
        {
            std::ostringstream e;
            e << fn << " does not have " << r+1 << " lines";
            throw exception(e);
        }
        return (*this)[r].size();
    }
    const std::string& operator()(unsigned i, unsigned j) const
    {
        std::ostringstream e;
        if(i>=(*this).size())
        {
            e << fn << " does not have " << i+1 << " lines";
            throw exception(e);
        }
        if(j>=(*this)[i].size())
        {
            e << "Line " << i+1
              << " in file " << fn
              << " does not have " << j+1
              << " columns (column " << c2excel(j) << ").";
            throw exception(e);
        }
        return (*this)[i][j];
    }

    unsigned getunsigned(unsigned i, unsigned j) const
    {
        unsigned res;
        std::string c = (*this)(i,j);
        try {
              res = stoul(c);
          } catch (...) {
              std::ostringstream e;
              e << "Cannot convert '" << c <<
                   "' to unsigned in file " << fn
                << ", " << c2excel(i,j);
              throw exception(e);
          }
        return res;
    }

    double getdouble(unsigned i, unsigned j) const
    {
        double res;
        std::string c = (*this)(i,j);
        try {
              res = stod(c);
          } catch (...) {
              std::ostringstream e;
              e << "Cannot convert '" << c
                << "' to double in file "
                << fn << ", " << c2excel(i,j);
              throw exception(e);
          }
        return res;
    }
    unsigned findinrow(unsigned r,const std::string& s)
    {
        for(unsigned i=0; i<this->c(r); i++)
        {
            if((*this)(r,i).find(s)==0)
                return i;
        }
        std::ostringstream os;
        os << "String '" << s
              << "' cannot be found in row " << r+1
              << " of file " << fn << std::endl;
        throw exception(os);
    }
    csv(const std::string& afn) :
        std::vector<std::vector<std::string>>(getvecs(afn)),
      fn(afn)
    {}
private:
    const std::string fn;
};

}; // namespace

#endif // CSV_HPP
