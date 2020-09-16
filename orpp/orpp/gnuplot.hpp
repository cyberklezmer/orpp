#ifndef ORPP_TOOLS_GNUPLOT_HPP
#define ORPP_TOOLS_GNUPLOT_HPP
 


/// A class encapsulting calls of gnuplot

namespace orpp
{

class gnuplot
{
/*        string name;
        string fn;
        ofstream sc;
public:
        ofstream& script() { return sc; }
        const string datfn() { return fn+".dat"; }
        const string scriptfn() { return fn+".sc"; }
        const string getfn() { return fn; }
        gnuplot(const string& aname) :
          name(aname),
          fn(hfd::interdir() + "/gnuplot/" + aname),
          sc(scriptfn().c_str())
      {
                 if(!sc) //fixme - should not be done in constructor
                 {
                         cout << "Error opening script " << scriptfn() << endl;
                         throw 1;
                 }
                 sc << "set terminal postscript landscape solid " << endl
                                << "set output '" << analysis::latexdir() + "/" + aname + ".eps" << "'" << endl;
          }
        void process()
        {
                sc.flush();

                string cmd = kgnuplotpath + " " + scriptfn();
                if(system(cmd.c_str()))
                {
                        cerr << "Error executing " << cmd << endl;
                        throw 1;
                }
        }
*/
};

} // namespace
#endif


