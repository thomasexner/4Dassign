void printhelp(string errormessage){
	cout << errormessage << endl;	
	cout << endl;
	cout << "Usage: assignment [-h][-p <file>][-n <file>][-l <file>][-m <int>][-t <int>][-s <float>]" << endl;
	cout << "    [-h]          : This helpscreen" << endl;
	cout << "    [-p <file>]   : Pdb-file (default: \"1D3Z.pdb\")" << endl;
	cout << "    [-n <file>]   : NOESY peaklist (default: \"NOESY.peaks\")" << endl;
	cout << "    [-l <file>]   : Log-file (default: \"logfile\")" << endl;
	cout << "    [-c <file>]   : Chemical shifts file (optional parameter: no default value)" << endl;
	cout << "    [-f <file>]   : Fix chemical shifts (default: no default value)" << endl;
	cout << "    [-s <float>]  : Safety distance added to restraints (default: 0.2)" << endl;
	cout << "    [-m <int>]    : Maximum number of Cycles (default: 100)" << endl;
	cout << "    [-t <int>]    : Number of threads used by openmp (default: 4)" << endl;
	cout << "    [-r <bool>]   : Call referencing routine, 0 = FALSE, 1 = TRUE (default: 0)" << endl;
	cout << "    [-d <bool>]   : Call identify doublepeaks routine, 0 = FALSE, 1 = TRUE (default: 0)" << endl;
	

	
	cout << endl;
	cout << endl;

	exit(1);
}

float string2float(string s){
	stringstream ss;
	float f = 0;
	
	ss << s;
	ss >> f;
	
	return f;
}

int string2int(string s){
	stringstream ss;
	int i = 0;
	
	ss << s;
	ss >> i;
	
	return i;
}
unsigned int string2uint(string s){
	stringstream ss;
	unsigned int i = 0;
	
	ss << s;
	ss >> i;
	
	return i;
}


vector<string> readfile(string file, bool die){
	vector<string> data;
	ifstream input;
	string line;
	
	input.open(file.c_str());
	if (input.good() == false){
		cerr << "ERROR! Cannot read: " << file << endl;
		input.close();
		input.clear();
		if (die == true) exit(2);
		if (die == false) return data;
	}
	while (1){
		getline(input,line);
		if (line.length() > 0) data.push_back(line);
		if (input.eof() == true) break;
	}
	input.close();
	input.clear();
	
	return data;
	
}

vector<atom> readatoms(string file){
	vector<atom> atoms;
	vector<string> data;
	
	data = readfile(file);
	
	for (unsigned int i = 0; i < data.size(); i++){
		if (data[i].length() < 4) continue;
		if (data[i].substr(0,4) != "ATOM") continue;
		
		atom peak;
		
		peak.name = data[i].substr(12,4);
		peak.res = data[i].substr(17,3);
		peak.resnr = string2uint(data[i].substr(23,3));
		peak.x = string2float(data[i].substr(31,7));
		peak.y = string2float(data[i].substr(39,7));
		peak.z = string2float(data[i].substr(47,7));
		peak.h = 0;
		peak.n = 0;
		peak.assigned = false;
		
		while (peak.name[0] == ' ') peak.name = peak.name.substr(1, (peak.name.length() - 1));
		while (peak.name[(peak.name.length() - 1)] == ' ') peak.name = peak.name.substr(0, (peak.name.length() - 1));

		atoms.push_back(peak);
	}
	
	// Remove all hydrogens that are not backbone NH
	for (unsigned int i = 0; i < atoms.size(); i++){
		if (atoms[i].name == "H") continue;

		atoms.erase(atoms.begin() + i);
		i--;
	}
	
	return atoms;
}


vector<atom> readfixed(string file){
	vector<atom> fix;
	vector<string> data;
	
	data = readfile(file);
	for (unsigned int i = 0; i < data.size(); i++){
		stringstream ss;
		float f1;
		float f2;
		unsigned int i1;
		unsigned int i2;
		string s1;
		atom peak;
		
		if (data[i].length() == 0) continue; 
		if (data[i][0] == '#') continue;
		ss << data[i];
		ss >> i1 >> i2 >> s1 >> f1 >> f2;
		// "1 21 ALA 7.355 212.55"

		peak.res = s1;
		peak.resnr = i2;
		peak.n = f1;
		peak.h = f2;
		
		fix.push_back(peak);
	}
	
	return fix;

}

vector<noe> readnoes(string file, bool ref){
	vector<noe> noes;
	vector<string> data;
	float reference = 0;
	
	data = readfile(file);

	// Use # for comments in NOESY file
	// and use entry with all zeros for reference: see NOESY file
	for (unsigned int i = 0; i < data.size(); i++){
		if (data[i].length() < 1) continue;
		if (data[i][0] == '#') continue;
		stringstream ss;
		float f1;
		float f2;
		float f3;
		float f4;
		int i1;
		int i2;
		int i3;
		ss << data[i];
		ss >> f1 >> f2 >> f3 >> f4 >> i1 >> i2 >> i3;
		if (abs(f1) < 0.1 && abs(f2) < 0.1 && abs(f3) < 0.1 && abs(f4) < 0.1) reference = float(i1);
		else {
			noe peak;
			peak.h1 = f1;
			peak.n1 = f2;
			peak.h2 = f3;
			peak.n2 = f4;
			peak.volume = i1;
			peak.doublepeak = false;
			peak.doublepeak2 = false;
			peak.dist = 0;
			peak.assigned = false;
			peak.nh2 = false;

			if (i3 == 1) peak.doublepeak = true;
			noes.push_back(peak);
		}
	}
	
	if (ref == true) return noes;
	
	// Generate distance restraints by using the reference:
	reference = reference * 4096;
	for (unsigned int i = 0; i < noes.size(); i++){
		noes[i].dist = 7.5;
		if (noes[i].volume < 1) continue;
		float distance = float(noes[i].volume);
		distance = reference / distance;
		distance = pow(distance, 0.16666666);
		noes[i].dist = distance;
		
	}
	
	return noes;
}


float referencing(vector<atom> & atoms, string noefile, string logfile, float safety, vector<atom> fix){
	unsigned int errlvl = 0;
	vector<noe> noes;
	float maxvol = 0;
	float reference = 0;
	
	noes = readnoes(noefile, true);
	for (unsigned int i = 0; i < noes.size(); i++){
		if (abs(noes[i].h1 - noes[i].h2) < 0.1 && abs(noes[i].n1 - noes[i].n2) < 0.1) continue;
		if (float(noes[i].volume) > maxvol) maxvol = float(noes[i].volume);
	}
	cout << "Referencing on peak volume: " << maxvol << endl;
	for (unsigned int i = 1; i < 50; i++){
		
		reference = float(i) * 0.1;
		
		reference = pow(reference, 6);
		reference = maxvol * reference; // this is A in V = A r^(-6)
		if (reference < 4096) continue;
		
		
		
		for (unsigned int j = 0; j < noes.size(); j++){
			noes[j].dist = 7.5;
			if (noes[j].volume < 1) continue;
			float distance = float(noes[j].volume);
			distance = reference / distance;
			distance = pow(distance, 0.16666666);
			if (distance > 7.5) continue;
			noes[j].dist = distance;
		}
		errlvl = optimize(atoms, noes, logfile, safety, fix);	
		// cout << "Tried: " << int(reference / 4096.0) << " => " << errlvl << endl;
		if (errlvl != 0) continue;
		cout << "Found valid restraint value: " << int((reference / 4096.0) + 0.5) << endl;
		return (reference / 4096.0);
		
		
	}
	
	return 0;
}

