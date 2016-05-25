// Place all used structures and classes here:
struct noe{
	float h1;
	float h2;
	float n1;
	float n2;
	float dist;
	float score;
	int volume;
	bool assigned;
	bool doublepeak;
	bool doublepeak2;
	bool nh2;
	vector<unsigned int> atoms;
	vector<unsigned int> atoms2;
	vector<unsigned int> crosspeaks;
};

struct atom{
	string name;
	string res;
	unsigned int resnr;
	float x;
	float y;
	float z;
	float n;
	float h;
	bool assigned;
	vector<unsigned int> noes;
};

struct statistic{
	float n;
	float h;
   vector<string> res;
   vector<unsigned int> resnr;
   vector<unsigned int> occurence;
};


// Place all routines here;
// prepare.cc:
void printhelp(string errormessage);
float string2float(string s);
int string2int(string s);
unsigned int string2uint(string s);
vector<string> readfile(string file, bool die = false);
vector<atom> readatoms(string file);
vector<atom> readfixed(string file);
vector<noe> readnoes(string file, bool ref = false);
float referencing(vector<atom> & atoms, string noefile, string logfile, float safety, vector<atom> fix);

// routines.cc:
unsigned int optimize(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, vector<atom> fix);
unsigned int artefacts(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, vector<atom> fix);
unsigned int examineassignment(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety);
unsigned int identifydp(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, vector<atom> fix);
bool checkassignment(vector<atom> atoms, vector<noe> noes, unsigned int atompeak, unsigned int noepeak, float safety);
void assignnh2(vector<atom> & atoms, vector<noe> & noes, float safety);
bool checkout(vector<vector<unsigned int> > assignment);
void getmatrix(vector<vector<unsigned int> > assignment, bool & check);
void assignall(vector<atom> & atoms, vector<noe> & noes, float safety, string logfile, vector<atom> fix);
void reduceassignment(vector<atom> & atoms, vector<noe> & noes, float safety);
bool checkdoublepeak(vector<atom> atoms, vector<noe> noes, unsigned int noe0, unsigned int noe1, unsigned int noe2, unsigned int atom00, unsigned int atom01, unsigned int atom10, unsigned int atom11, float safety);
bool checkreduced(vector<atom> atoms, vector<noe> noes, float safety);
void speedup(vector<noe> & noes, vector<atom> atoms);
void fix_atoms(vector<atom> & atoms, vector<noe> & noes, vector<atom> fix);


// results.cc
void printassignment(vector<atom> atoms, vector<noe> noes, string logfile);
void printnoes(vector<atom> atoms, vector<noe> noes, string logfile);
void getstatistics(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, string csfile);
void reduceassis(vector<noe> noes, vector<vector<noe> > & allassis, vector<atom> & atoms, float & safety);
void writestatistics(vector<atom> & atoms, vector<noe> & noes, string logfile, vector<vector<noe> > & allassis);
void scoreassis(vector<vector<noe> > & allassis, vector<atom> atoms, vector<noe> noes, float safety, string logfile);
void shiftassis(vector<vector<noe> > & allassis, vector<atom> atoms, vector<noe> noes, string csfile, string logfile);
