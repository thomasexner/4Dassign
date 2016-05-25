// Never use any other namespace than std! Leave this always as first lines!
using namespace std;   

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>

#include "header.h"
#include "prepare.cc"
#include "routines.cc"
#include "results.cc"



int main (int argc, char* argv[]){
	vector<atom> atoms;
	vector<atom> fix;
	vector<noe> noes;
	unsigned int errorlevel = 0;
	float safety = 0.2;
	ofstream output;
	bool reference = false;
	bool identdp = false;

	// Initialize paramters with default values:
	string pdbfile = "1D3Z.pdb";
	string noefile = "NOESY.peaks";
	string logfile = "logfile";
	string csfile;
	string fixed_file;

	// Initialize random generator for all subroutines
	srand(time(NULL));

	// read parameters, if a parameter is missing print helpscreen and exit
	for (int i = 0; i < argc; i++){
		string s = argv[i];
		
		if (s == "-h" || s == "--help") printhelp("");
		
		if (s == "-p" && argc > i) pdbfile = argv[(i + 1)];
		if (s == "-p" && argc == i) printhelp("Missing expression after -p parameter");

		if (s == "-n" && argc > i) noefile = argv[(i + 1)];
		if (s == "-n" && argc == i) printhelp("Missing expression after -n parameter");
		
		if (s == "-l" && argc > i) logfile = argv[(i + 1)];
		if (s == "-l" && argc == i) printhelp("Missing expression after -l parameter");

		if (s == "-c" && argc > i) csfile = argv[(i + 1)];
		if (s == "-c" && argc == i) printhelp("Missing expression after -c parameter");

		if (s == "-s" && argc > i) safety = string2float(argv[(i + 1)]);
		if (s == "-s" && argc == i) printhelp("Missing expression after -s parameter");

		if (s == "-r" && argc > i) if (string2int(argv[(i + 1)]) != 0) reference = true;
		if (s == "-r" && argc == i) reference = true;

		if (s == "-d" && argc > i) if (string2int(argv[(i + 1)]) != 0) identdp = true;
		if (s == "-d" && argc == i) identdp = true;
		
		if (s == "-f" && argc > i) fixed_file = argv[(i + 1)];
		if (s == "-f" && argc == i) printhelp("Missing expression after -f parameter");
		
		
		
		

	}
	
	
	
	if (reference == true){
		cout << "Start referencing routine" << endl;
		float reference = 0.0;

		atoms = readatoms(pdbfile);
		reference = referencing(atoms, noefile, logfile, safety, fix);		

		output.open(logfile.c_str());
		output << "Reference value found: " << int(reference + 0.5) << endl;
		output.close();
		output.clear();
		return 0;
		
	}
	// Read pdbfile and remove all heavy atoms and all non-amin hydrogens 
	atoms = readatoms(pdbfile);
	noes = readnoes(noefile);
	if (fixed_file.length() > 0) fix = readfixed(fixed_file);

	// Create an empty logfile
	output.open(logfile.c_str());
	output.close();
	output.clear();
	if (identdp == true){
		errorlevel = identifydp(atoms, noes, logfile, safety, fix);	
		return errorlevel;
	}
	while(1){
		errorlevel = optimize(atoms, noes, logfile, safety, fix);	
		if (errorlevel == 0){
			cout << "Start Referencing routine" << endl;
			getstatistics(atoms, noes, logfile, safety, csfile);	
			return 0;
			break;
		}
		else {
			cout << "Start artefacts routine" << endl;
			vector<noe> noecopy = noes;
			errorlevel = artefacts(atoms, noes, logfile, safety, fix);
			if (errorlevel != 0){
				cerr << "Unable to identify artefacts: Errorcode: " << errorlevel << endl;
				cerr << "Exiting." << endl;
				return errorlevel;
			}
			if (errorlevel == 0) getstatistics(atoms, noes, logfile, safety, csfile);	
			return 0;
			break;
		}
	}
	
	return 0;
	
}
