#ifndef ORPP_GNUPLOT_HPP
#define ORPP_GNUPLOT_HPP
 
#include "orpp.hpp"

/// A class encapsulting calls of gnuplot

namespace orpp
{

class gnuplot
{
    std::string name;
    std::ofstream sc;
    std::ofstream dat;
    static constexpr const char* klabel = "gnuplot";
    static std::string scriptfn(const std::string& an)
       { return sys::tmpfolder() + an + ".plt"; }
    static std::string datfn(const std::string& an)
       { return sys::tmpfolder() + an + ".dat"; }
public:
    const std::string datfn() { return datfn(name); }
    std::ofstream& script() { return sc; }
    std::ofstream& datfile() { return dat; }
    gnuplot(const std::string& aname) :
      name(aname), sc(scriptfn(aname)), dat(datfn(aname))
      {
         sc << "set terminal postscript eps color " << std::endl;
         sc << "set output '"
            << sys::outputfolder() + aname + ".eps" << "'" << std::endl;
      }
    void process()
    {
        sc.flush();
        if(!sc)
            throw "Error was writing to script " + scriptfn(name);
        dat.flush();
        if(!dat)
            throw "Error was writing to data file " + datfn(name);
        std::string cmd = sys::get(klabel).second + klabel + " " + scriptfn(name);
        if(system(cmd.c_str()))
            throw "Error executing " + cmd;
    }
};


} // namespace
#endif


