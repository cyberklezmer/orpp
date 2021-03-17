
#ifndef ORPP_TOOLS_GRETL_HPP
#define ORPP_TOOLS_GRETL_HPP
 


/// A class encapsulting calls of GRETL

class gretl
{
	string dir;
	string tmpfn()
	{
		return dir+"/_gretl_tmp";
	}
	string scriptresult;
public:
	gretl(const string& adir): dir(adir) {}
	string scdir()
	{
		return dir;
	}
	void runscript(const string& scfn)
	{
		using namespace std;
		string cmd = "gretlcli -b \"" + scfn + "\" > " + tmpfn();
		cout << cmd << endl;
		if(system(cmd.c_str()))
		{
			cerr << "Failed to run " << cmd << endl;
			throw 1;
		}
		ifstream p(tmpfn().c_str());
		if(!p)
		{
			cerr << "Error opening " << tmpfn() << endl;
			throw 1;
		}

		getline(p,scriptresult,'\0');
	}
	string getscriptoutput()
	{
		ifstream s(tmpfn().c_str());
		if(!s)
		{
                        cerr << "Error opening " << tmpfn().c_str() << endl;
			throw 1;
		}
		string str;
		getline(s,str);
		return str;
	}
	double findpar(const string& parname) // assumes that "? parname" appears in the output
	{
		string label = "? " + parname;
		unsigned int pos;
		if((pos = scriptresult.find(label))==string::npos)
		{
		    cerr << "Cannot find parameter " << parname << " in gretl script output " << endl;
			throw 1;
		}
		while(scriptresult[pos++] != '\n')
			;
		return atof(scriptresult.c_str()+pos);
	}
};

#endif


