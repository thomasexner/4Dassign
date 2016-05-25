void printassignment(vector<atom> atoms, vector<noe> noes, string logfile){
	ofstream output;
	vector<atom> exps;
	
	// Read experimental shifts and store them in exps vector
	vector<string> data = readfile("exp.nmr");
	for (unsigned int i = 0; i < data.size(); i++){
		stringstream ss;
		string s1;
		
		string s2;
		float f1;
		int i1;
		int i2;
		atom peak;
		bool accept = false;
		
		ss << data[i];
		ss >> i1 >> i2 >> s1 >> s2 >> f1;
		
		if (s2 == "H") accept = true;
		if (s1 == "ASN" && s2 == "HD21") accept = true;
		if (s1 == "ASN" && s2 == "HD22") accept = true;
		if (s1 == "GLN" && s2 == "HE21") accept = true;
		if (s1 == "GLN" && s2 == "HE22") accept = true;
		if (s1 == "TRP" && s2 == "HE1") accept = true;
		if (s1 == "ARG" && s2 == "HE") accept = true;
		
		if (accept == false) continue;
		
		peak.res = s1;
		peak.resnr = i2;
		peak.name = s2;
		peak.h = f1;
		peak.n = 0;
		peak.x = 0;
		peak.y = 0;
		peak.z = 0;
		peak.assigned = true;
		exps.push_back(peak);
	}

	for (unsigned int i = 0; i < data.size(); i++){
		stringstream ss;
		string s1;
		string s2;
		float f1;
		unsigned int i1;
		unsigned int i2;
		atom peak;
		
		ss << data[i];
		ss >> i1 >> i2 >> s1 >> s2 >> f1;
		
		if (s2 == "N") for (unsigned int j = 0; j < exps.size(); j++) if (exps[j].resnr == i2 && exps[j].name == "H") exps[j].n = f1;
		if (s1 == "ASN" && s2 == "ND2") for (unsigned int j = 0; j < exps.size(); j++) if (exps[j].resnr == i2) if (exps[j].name == "HD21" || exps[j].name == "HD22") exps[j].n = f1;
		if (s1 == "GLN" && s2 == "NE2") for (unsigned int j = 0; j < exps.size(); j++) if (exps[j].resnr == i2) if (exps[j].name == "HE21" || exps[j].name == "HE22") exps[j].n = f1;
		if (s1 == "TRP" && s2 == "NE1") for (unsigned int j = 0; j < exps.size(); j++) if (exps[j].resnr == i2 && exps[j].name == "HE1") exps[j].n = f1;
		if (s1 == "ARG" && s2 == "NE") for (unsigned int j = 0; j < exps.size(); j++) if (exps[j].resnr == i2 && exps[j].name == "HE") exps[j].n = f1 + 32;
		
	}
	data.clear();	

	
	output.open(logfile.c_str(), ios::app);
	output << "Current results for atoms: " << endl;
	
	for (unsigned int i = 0; i < atoms.size(); i++){
		if (atoms[i].resnr < 10) output << " ";
		output << atoms[i].resnr << " " << atoms[i].res << " ";
		if (atoms[i].name.length() < 4) output << " ";
		if (atoms[i].name.length() < 3) output << " ";
		if (atoms[i].name.length() < 2) output << " ";
		output << atoms[i].name << " => ";
		if (atoms[i].assigned == false){
			output << "n.a. ";
			for (unsigned int j = 0; j < atoms[i].noes.size(); j++){
				if (j > 0) output << "+ ";
				output << fixed << setprecision(2) << noes[(atoms[i].noes[j])].h1 << " ";
				output << fixed << setprecision(1) << noes[(atoms[i].noes[j])].n1 << " ";
			}
		}
		if (atoms[i].assigned == true) output << fixed << setprecision(2) << atoms[i].h << " " << fixed << setprecision(1) << atoms[i].n << " ";
		output << "(";
		bool found = false;
		for (unsigned int j = 0; j < exps.size(); j++) if (exps[j].resnr == atoms[i].resnr && atoms[i].name == exps[j].name){
			found = true;
			output << fixed << setprecision(2) << exps[j].h << " " << fixed << setprecision(1) << exps[j].n <<  ") ";
		}
		if (found == false) output << "n.a.)       ";
		output << endl;
	}
	
	output.close();
	output.clear();

	return;
}

void printnoes(vector<atom> atoms, vector<noe> noes, string logfile){
	ofstream output;
	
	output.open(logfile.c_str(),ios::app);
	output << "Overview over all NOEs:" << endl;
	output << "shifts (number of crosspeaks) => remaining assignment possibilities" << endl;
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h2 - noes[i].h1) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1){
		output << fixed << setprecision(2) << noes[i].n1 << " " << noes[i].h1 << " (" << noes[i].crosspeaks.size() << ") => ";
		for (unsigned int j = 0; j < noes[i].atoms.size(); j++) output << atoms[(noes[i].atoms[j])].resnr << " " << atoms[(noes[i].atoms[j])].res << " " << atoms[(noes[i].atoms[j])].name << " ";
		output << endl;
	}
	
	output.close();
	output.clear();
	
	return;
}

void getstatistics(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, string csfile){
	// Write the final (still ambigous) assignment to the logfile
	// printassignment(atoms, noes, logfile);
	printnoes(atoms, noes, logfile);
	
	exit(1);
	
	vector<vector<noe> > allassis;
	vector<float> scores;
	
	// **************************
	// remove double peaks here cannot handle them properly

	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].doublepeak == true){
		noe peak = noes[i];
		for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1){
			noes.erase(noes.begin() + j);
			if (j < i) i--;
			j--;
		}
	
		
		for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[i].h1 - noes[j].h2) < 0.01 && abs(noes[i].n1 - noes[j].n2) < 0.1){
			noes.erase(noes.begin() + j);
			if (j < i) i--;
			j--;
		}
		
				
		
		
		noe peak1 = peak;
		noe peak2 = peak;

		noes.erase(noes.begin() + i);
		i--;
		
		peak1.doublepeak = false;
		peak2.doublepeak = false;
		peak1.atoms.clear();
		peak2.atoms.clear();
		peak1.crosspeaks.clear();
		peak2.crosspeaks.clear();
		peak1.assigned = true;
		peak2.assigned = true;
		peak1.atoms.push_back(peak.atoms[0]);
		peak2.atoms.push_back(peak.atoms[1]);
		peak2.n1 = peak2.n1 + 0.11;
		peak2.n2 = peak2.n2 + 0.11;
		peak2.h1 = peak2.h1 + 0.1;
		peak2.h2 = peak2.h2 + 0.1;
		noes.push_back(peak1);
		noes.push_back(peak2);

	
	}
	
	// Generate all possible assignments out of this assignment matrix
	// cout << "Generating assignment possibilities" << endl;
	
	// Build new crosspeaks vector
	for (unsigned int i = 0; i < noes.size(); i++) noes[i].crosspeaks.clear();
	for (unsigned int i = 0; i < noes.size(); i++) for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1) noes[i].crosspeaks.push_back(j);
		
	// reduceassis(noes, allassis, atoms, safety);
	
   if (allassis.size() > 0){
      cout << "Safe assignment to 'tmp'" << endl;
      ofstream output;
     
      output.open("tmp",ios::app);
      for (unsigned int i = 0; i < allassis.size(); i++){ 
			for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1) {
   	  	 	if (allassis[i][j].atoms.size() != 1) output << "(";
   	   	for (unsigned int k = 0; k < allassis[i][j].atoms.size(); k++) output <<  allassis[i][j].atoms[k] << " ";
   	   	if (allassis[i][j].atoms.size() != 1) output << ")";
				output << " ";
			}
   	   output << endl;
      }
      output.close();
      output.clear();
   }
	
	// exit(0);
	cout << "Read assignment possibilities and call routine" << endl;
	ofstream output;
	
	while(1){
		ifstream input;
		string line;
		
		input.open("tmp");
		if (input.good() == false){
			cerr << "ERROR! Cannot read 'tmp'" << endl;
			input.close();
			input.clear();
			exit(2);
		}
		while (1){
			getline(input,line);
			if (line.length() == 0) break;
			vector<noe> peak = noes;
			for (unsigned int i = 0; i < peak.size(); i++) if (abs(peak[i].h1 - peak[i].h2) < 0.01 && abs(peak[i].n1 - peak[i].n2) < 0.1){
				if (line.length() < 1) break;
				while(line[0] == ' ' && line.length() > 1) line = line.substr(1, (line.length() - 1));
				if (line.length() < 1) break;
				for (unsigned int j = 0; j < line.length(); j++) if (line[j] == ' '){
					string nr = line.substr(0,j);
					// cout << " " << nr << " -> " << string2float(nr) << " ";
					peak[i].atoms.clear();
					peak[i].atoms.push_back(string2float(nr));
					while(line[0] != ' ' && line.length() > 1) line = line.substr(1, (line.length() - 1));
					while(line[0] == ' ' && line.length() > 1) line = line.substr(1, (line.length() - 1));
				break;
				}
				if (line.length() < 1) break;
			}
						
			
			allassis.push_back(peak);
			if (allassis.size() % 100 == 0) cout << ".";
			if (allassis.size() >= 1000){
				cout << endl;
				cout << "Remove higher scored entries (" << allassis.size() << ") " << endl;
				vector<long double> scores;
				
				for (unsigned int i = 0; i < allassis.size(); i++) for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
					cout << fixed << setprecision(2) << allassis[i][j].h1 << " " << allassis[i][j].n1 << " - " << allassis[i][j].atoms[0] << atoms[(allassis[i][j].atoms[0])].resnr << " ";
					cout << atoms[(allassis[i][j].atoms[0])].res << " => (" << allassis[i][j].crosspeaks.size() << ") ";
					for (unsigned int k = 0; k < allassis[i][j].crosspeaks.size(); k++) for (unsigned int m = 0; m < allassis[i].size(); m++) if (abs(allassis[i][m].h1 - allassis[i][m].h2) < 0.01 && abs(allassis[i][m].n1 - allassis[i][m].n2) < 0.1){
						if (m != j) if (abs(allassis[i][m].h1 - allassis[i][(allassis[i][j].crosspeaks[k])].h2) < 0.01 && abs(allassis[i][m].n1 - allassis[i][(allassis[i][j].crosspeaks[k])].n2) < 0.1)
						cout << atoms[(allassis[i][m].atoms[0])].resnr << " " << atoms[(allassis[i][m].atoms[0])].res << " ";
					}
					cout << endl;
				}
				
				for (unsigned int i = 0; i < allassis.size(); i++){
					long double dist = 0;
					for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1) for (unsigned int k = 0; k < allassis[i][j].crosspeaks.size(); k++){
						unsigned int diag = 0;	
						long double score = 0;
							
						for (unsigned int m = 0; m < allassis[i].size(); m++) if (abs(allassis[i][m].h1 - allassis[i][m].h2) < 0.01 && abs(allassis[i][m].n1 - allassis[i][m].n2) < 0.1){
							if (abs(allassis[i][(allassis[i][j].crosspeaks[k])].h2 - allassis[i][m].h1) < 0.01 && abs(allassis[i][(allassis[i][j].crosspeaks[k])].n2 - allassis[i][m].n1) < 0.1) diag = m;
						}
						score = score + ((atoms[(allassis[i][j].atoms[0])].x - atoms[(allassis[i][diag].atoms[0])].x) * (atoms[(allassis[i][j].atoms[0])].x - atoms[(allassis[i][diag].atoms[0])].x));
						score = score + ((atoms[(allassis[i][j].atoms[0])].y - atoms[(allassis[i][diag].atoms[0])].y) * (atoms[(allassis[i][j].atoms[0])].y - atoms[(allassis[i][diag].atoms[0])].y));
						score = score + ((atoms[(allassis[i][j].atoms[0])].z - atoms[(allassis[i][diag].atoms[0])].z) * (atoms[(allassis[i][j].atoms[0])].z - atoms[(allassis[i][diag].atoms[0])].z));
			
						score = sqrt(score);
						score = abs(score -  noes[(allassis[i][j].crosspeaks[k])].dist);
						dist = dist + (score * score);
					}
					// cout << dist << endl;
					scores.push_back(dist);	
				}
				vector<long double> ranks = scores;
				long double cutoff;

				sort(ranks.begin(), ranks.end());
				
				cutoff = ranks[10];				
				cout << cutoff << endl;
				
				if (scores.size() != allassis.size()){
					cerr << "ERROR in calculation of scores" << endl;
					exit(1);
				}
				output.open("out.tmp",ios::app);
				for (unsigned int i = 0; i < scores.size(); i++) if (scores[i] > cutoff){
					// cout << "Remove entry" << endl;

					output << fixed << setprecision(1) << scores[i] << " ";
					for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
						string res = atoms[(allassis[i][j].atoms[0])].res;
					 	if (res == "GLY") output << "G";
					 	if (res == "ALA") output << "A";
						if (res == "SER") output << "S";
						if (res == "THR") output << "T";
					 	if (res == "GLN") output << "Q";
						if (res == "GLU") output << "E";
						if (res == "ASN") output << "N";
						if (res == "ASP") output << "D";
					 	if (res == "CYS") output << "C";
					 	if (res == "MET") output << "M";
					 	if (res == "VAL") output << "V";
					 	if (res == "ILE") output << "I";
					 	if (res == "LEU") output << "L";
					 	if (res == "ARG") output << "R";
					 	if (res == "LYS") output << "K";
					  	if (res == "HIS") output << "H";
					  	if (res == "TYR") output << "Y";
					 	if (res == "TRP") output << "W";
					 	if (res == "PHE") output << "F";
					 	output << atoms[(allassis[i][j].atoms[0])].resnr << " ";
					}
					output << endl;

					
					scores.erase(scores.begin() + i);
					allassis.erase(allassis.begin() + i);
					i--;
				}
				output.close();
				output.clear();
				
				cout << "Done: (" << allassis.size() << ")" << endl;
			}
			if (input.eof() == true) break;
		}
		input.close();
		input.clear();

		break;
	}
	
	output.open("out.tmp",ios::app);
	for (unsigned int i = 0; i < scores.size(); i++){
		// cout << "Remove entry" << endl;

		output << fixed << setprecision(1) << scores[i] << " ";
		for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
			string res = atoms[(allassis[i][j].atoms[0])].res;
		 	if (res == "GLY") output << "G";
		 	if (res == "ALA") output << "A";
			if (res == "SER") output << "S";
			if (res == "THR") output << "T";
		 	if (res == "GLN") output << "Q";
			if (res == "GLU") output << "E";
			if (res == "ASN") output << "N";
			if (res == "ASP") output << "D";
		 	if (res == "CYS") output << "C";
		 	if (res == "MET") output << "M";
		 	if (res == "VAL") output << "V";
		 	if (res == "ILE") output << "I";
		 	if (res == "LEU") output << "L";
		 	if (res == "ARG") output << "R";
		 	if (res == "LYS") output << "K";
		  	if (res == "HIS") output << "H";
		  	if (res == "TYR") output << "Y";
		 	if (res == "TRP") output << "W";
		 	if (res == "PHE") output << "F";
		 	output << atoms[(allassis[i][j].atoms[0])].resnr << " ";
		}
		output << endl;
	}
	output.close();
	output.clear();
	
	
	exit(0);
	
	// return;
	
	// Give the remaining assignment a score, depending on how well the fulfill the restraints, i.e. the RMSD of the restraints not just an upper distance restraint.
	// cout << "There are " << allassis.size() << " assignment possibilities left." << endl;
	// writestatistics(atoms, noes, logfile, allassis);
	// cout << "Score remaining assignments by their distance restraints." << endl;
	// scoreassis(allassis, atoms, noes, safety, logfile);
	if (csfile.length() > 0){
		vector<string> data;
		vector<atom> shifts;
	
		data = readfile(csfile);
		for (unsigned int i = 0; i < data.size(); i++){
			stringstream ss;
			unsigned int i1;
			unsigned int i2;
			string s1;
			string s2;
			float f1;
			atom peak;
		
			ss << data[i];
			ss >> i1 >> i2 >> s1 >> s2 >> f1;
			if (s2 != "H") continue;
			peak.resnr = i2;
			peak.res = s1;
			peak.h = f1;
			shifts.push_back(peak);
		}


		for (unsigned int i = 0; i < data.size(); i++){
			stringstream ss;
			unsigned int i1;
			unsigned int i2;
			string s1;
			string s2;
			float f1;
		
			ss << data[i];
			ss >> i1 >> i2 >> s1 >> s2 >> f1;
			if (s2 != "N") continue;
			for (unsigned int j = 0; j < shifts.size(); j++) if (shifts[j].res == s1 && shifts[j].resnr == i2) shifts[j].n = f1;		
		}		
		
		ifstream input;
		string line;
		input.open("tmp");
		if (input.good() == false){
			cerr << "ERROR! Cannot read 'tmp'" << endl;
			input.close();
			input.clear();
			exit(2);
		}
		while (1){
			getline(input,line);
			if (line.length() == 0) break;
			vector<noe> peak = noes;
			for (unsigned int i = 0; i < peak.size(); i++) if (abs(peak[i].h1 - peak[i].h2) < 0.01 && abs(peak[i].n1 - peak[i].n2) < 0.1){
				if (line.length() < 1) break;
				while(line[0] == ' ' && line.length() > 1) line = line.substr(1, (line.length() - 1));
				if (line.length() < 1) break;
				for (unsigned int j = 0; j < line.length(); j++) if (line[j] == ' '){
					string nr = line.substr(0,j);
					// cout << " " << nr << " -> " << string2float(nr) << " ";
					peak[i].atoms.clear();
					peak[i].atoms.push_back(string2float(nr));
					while(line[0] != ' ' && line.length() > 1) line = line.substr(1, (line.length() - 1));
					while(line[0] == ' ' && line.length() > 1) line = line.substr(1, (line.length() - 1));
				break;
				}
				if (line.length() < 1) break;
			}
					
			
			allassis.push_back(peak);
			if (allassis.size() % 100 == 0) cout << ".";
			if (allassis.size() >= 1000){
				cout << endl;
				cout << "Remove higher scored entries (" << allassis.size() << ") " << endl;
				vector<float> scores;
				
				for (unsigned int i = 0; i < allassis.size(); i++){
					float score = 0;
					for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
						float value = 0;
						for (unsigned int k = 0; k < shifts.size(); k++) if (shifts[k].resnr == atoms[(allassis[i][j].atoms[0])].resnr){
							// cout << allassis[i][j].n1 << " - " << shifts[k].n << endl;
							value = abs(allassis[i][j].n1 - shifts[k].n);
							value = value * 10;
							value = value + abs(allassis[i][j].h1 - shifts[k].h);
						}
						score = score + value;
					}
					scores.push_back(score);
				}
				vector<float> ranks = scores;
				float cutoff;

				sort(ranks.begin(), ranks.end());
				
				cutoff = ranks[10];				
				cout << cutoff << endl;
				if (scores.size() != allassis.size()){
					cerr << "ERROR in calculation of scores" << endl;
					exit(1);
				}
				ofstream output;
				output.open("out.tmp",ios::app);
				
				for (unsigned int i = 0; i < scores.size(); i++) if (scores[i] > cutoff){
					// cout << "Remove entry" << endl;
					
					output << fixed << setprecision(1) << scores[i] << " ";
					for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
						string res = atoms[(allassis[i][j].atoms[0])].res;
					 	if (res == "GLY") output << "G";
					 	if (res == "ALA") output << "A";
						if (res == "SER") output << "S";
						if (res == "THR") output << "T";
					 	if (res == "GLN") output << "Q";
						if (res == "GLU") output << "E";
						if (res == "ASN") output << "N";
						if (res == "ASP") output << "D";
					 	if (res == "CYS") output << "C";
					 	if (res == "MET") output << "M";
					 	if (res == "VAL") output << "V";
					 	if (res == "ILE") output << "I";
					 	if (res == "LEU") output << "L";
					 	if (res == "ARG") output << "R";
					 	if (res == "LYS") output << "K";
					  	if (res == "HIS") output << "H";
					  	if (res == "TYR") output << "Y";
					 	if (res == "TRP") output << "W";
					 	if (res == "PHE") output << "F";
					 	output << atoms[(allassis[i][j].atoms[0])].resnr << " ";
					}
					output << endl;


					scores.erase(scores.begin() + i);
					allassis.erase(allassis.begin() + i);
					i--;
				}
				output.close();
				output.clear();
				
				cout << "Done: (" << allassis.size() << ")" << endl;
			}
			if (input.eof() == true) break;
		}
		input.close();
		input.clear();
		shiftassis(allassis, atoms, noes, csfile, logfile);
	}
	return;
}

void writestatistics(vector<atom> & atoms, vector<noe> & noes, string logfile, vector<vector<noe> > & allassis){
	
	vector<statistic> statistics;
	ofstream output;
	
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1){
		statistic peak;
		peak.n = noes[i].n1;
		peak.h = noes[i].h1;
		statistics.push_back(peak);
		if (noes[i].doublepeak == true){
			peak.n = peak.n + 0.1;
			statistics.push_back(peak);
		} 
	}
	for (unsigned int i = 0; i < statistics.size(); i++) for (unsigned int j = 0; j < atoms.size(); j++){
		statistics[i].res.push_back(atoms[j].res);
		statistics[i].resnr.push_back(atoms[j].resnr);
		statistics[i].occurence.push_back(0);
	}
	
	for (unsigned int i = 0; i < statistics.size(); i++) for (unsigned int j = 0; j < allassis.size(); j++) for (unsigned int k = 0; k < allassis[j].size(); k++){ 
		if (abs(allassis[j][k].h1 - allassis[j][k].h2) < 0.01 && abs(allassis[j][k].n1 - allassis[j][k].n2) < 0.1) if (abs(allassis[j][k].h1 - statistics[i].h) < 0.01 && abs(allassis[j][k].n1 - statistics[i].n) < 0.01) if (allassis[j][k].atoms.size() == 1){
			for (unsigned int m = 0; m < statistics[i].resnr.size(); m++) if (statistics[i].resnr[m] == atoms[(allassis[j][k].atoms[0])].resnr) statistics[i].occurence[m] = statistics[i].occurence[m] + 1;
		}
	}
	
	
	output.open(logfile.c_str(),ios::app);
	output << "Statistics of remaining assignments:" << endl;
	for (unsigned int i = 0; i < statistics.size(); i++){	
		output << fixed << setprecision(2) << statistics[i].n << " " << statistics[i].h << " => ";
		for (unsigned int j = 0; j < statistics[i].occurence.size(); j++) if (statistics[i].occurence[j] > 0) output << statistics[i].resnr[j] << " " << statistics[i].res[j] << " (" << statistics[i].occurence[j] << ") ";
		output << endl;
		
	}
	output.close();
	output.clear();
}


void reduceassis(vector<noe> noes, vector<vector<noe> > & allassis, vector<atom> & atoms, float & safety){
	// first check if assignment is already reduced, i.e. unambigouasly (but not guaranteed correct)
	// if so validate assignment and store or reject assignment
	bool check = true;
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h2 - noes[i].h1) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.01) if (noes[i].atoms.size() > 1) check = false;
	if (check == true){
		vector<noe> peak;
		if (checkreduced(atoms, noes, safety) == false) return;
		
		check = true;
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h2 - noes[i].h1) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.01){
			for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[j].h2 - noes[j].h1) < 0.01 && abs(noes[j].n1 - noes[j].n2) < 0.01) if (noes[i].atoms[0] == noes[j].atoms[0]) check = false;
		}
		if (check == false) return;
		for (unsigned int i = 0; i < noes.size(); i++) peak.push_back(noes[i]);
		for (unsigned int i = 0; i < peak.size(); i++) if (abs(peak[i].h2 - peak[i].h1) < 0.01 && abs(peak[i].n1 - peak[i].n2) < 0.01) peak[i].assigned = true;
		if (checkreduced(atoms, peak, safety) == false) return;
		allassis.push_back(peak);
		if (allassis.size() % 10000 == 0){
			cout << "Safe assignment to 'tmp'" << endl;
			ofstream output;
			
      	output.open("tmp",ios::app);
   	   for (unsigned int i = 0; i < allassis.size(); i++){ 
				for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
   	  	 		if (allassis[i][j].atoms.size() != 1) output << "(";
   	   		for (unsigned int k = 0; k < allassis[i][j].atoms.size(); k++) output <<  allassis[i][j].atoms[k] << " ";
   		   	if (allassis[i][j].atoms.size() != 1) output << ")";
					output << " ";
				}
   	   	output << endl;
      	}
   	   output.close();
	      output.clear();
			allassis.clear();
			
		}
		return;
	}
	
	// If not yet reduced, go iteratively through noe-vector and find first ambigous entry
	// Generate as many sub-matrices as entries. Each submatrix contains one of the ambigous entries of the parent matrix
	// Call this routine with each submatrix and exit afterwards: All other ambigous entries are examined iterativly by the following calls of this routine with the submatrices
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h2 - noes[i].h1) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.01) if (noes[i].atoms.size() > 1){
		for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
			vector<noe> newnoes = noes;
			vector<atom> newatoms = atoms;
			newnoes[i].atoms.clear();
			newnoes[i].atoms.push_back(noes[i].atoms[j]);
			newnoes[i].assigned = true;
			for (unsigned int k = 0; k < noes.size(); k++) noes[k].crosspeaks.clear();
			for (unsigned int k = 0; k < noes.size(); k++) for (unsigned int m = 0; m < noes.size(); m++) if (m != k) if (abs(noes[m].h1 - noes[k].h1) < 0.01 && abs(noes[m].n1 - noes[k].n1) < 0.1) noes[k].crosspeaks.push_back(m);

			reduceassignment(newatoms, newnoes, safety);
			if (checkreduced(atoms, newnoes, safety) == false) continue;
			reduceassis(newnoes, allassis, atoms, safety);
		}
		return;
	}
	

	return;
}

void scoreassis(vector<vector<noe> > & allassis, vector<atom> atoms, vector<noe> noes, float safety, string logfile){
	vector<long double> scores;
	ofstream output;
	
	for (unsigned int i = 0; i < allassis.size(); i++){
		long double dist = 0;
		for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1) for (unsigned int k = 0; k < allassis[i][j].crosspeaks.size(); k++){
			unsigned int diag = 0;	
			long double score = 0;
						
			for (unsigned int m = 0; m < allassis[i].size(); m++) if (abs(allassis[i][m].h1 - allassis[i][m].h2) < 0.01 && abs(allassis[i][m].n1 - allassis[i][m].n2) < 0.1){
				if (abs(noes[(allassis[i][j].crosspeaks[k])].h2 - allassis[i][m].h1) < 0.01 && abs(noes[(allassis[i][j].crosspeaks[k])].n2 - allassis[i][m].n1) < 0.1) diag = m;
			}
			score = score + ((atoms[(allassis[i][j].atoms[0])].x - atoms[(allassis[i][diag].atoms[0])].x) * (atoms[(allassis[i][j].atoms[0])].x - atoms[(allassis[i][diag].atoms[0])].x));
			score = score + ((atoms[(allassis[i][j].atoms[0])].y - atoms[(allassis[i][diag].atoms[0])].y) * (atoms[(allassis[i][j].atoms[0])].y - atoms[(allassis[i][diag].atoms[0])].y));
			score = score + ((atoms[(allassis[i][j].atoms[0])].z - atoms[(allassis[i][diag].atoms[0])].z) * (atoms[(allassis[i][j].atoms[0])].z - atoms[(allassis[i][diag].atoms[0])].z));
			
			score = sqrt(score);
			score = abs(score -  noes[(allassis[i][j].crosspeaks[k])].dist);
			dist = dist + (score * score);
		}
		scores.push_back(dist);	
	}
	
	output.open(logfile.c_str(),ios::app);
	output << "Overview over ranked remaining assignments. Score is the sum over the absolute violations (in AA) of the restraints" << endl;
	for (unsigned int i = 0; i < scores.size(); i++){	
		float score = -1;
		unsigned int best = 0;
		for (unsigned int j = 0; j < scores.size(); j++) if (scores[j] > 0) if (score < 0 || scores[j] < score){
			score = scores[j];
			best = j;
		}
		output << fixed << setprecision(1) << score << " ";
		for (unsigned int j = 0; j < allassis[best].size(); j++) if (abs(allassis[best][j].h1 - allassis[best][j].h2) < 0.01 && abs(allassis[best][j].n1 - allassis[best][j].n2) < 0.1){
			 string res = atoms[(allassis[best][j].atoms[0])].res;
			 if (res == "GLY") output << "G";
			 if (res == "ALA") output << "A";
			 if (res == "SER") output << "S";
			 if (res == "THR") output << "T";
			 if (res == "GLN") output << "Q";
			 if (res == "GLU") output << "E";
			 if (res == "ASN") output << "N";
			 if (res == "ASP") output << "D";
			 if (res == "CYS") output << "C";
			 if (res == "MET") output << "M";
			 if (res == "VAL") output << "V";
			 if (res == "ILE") output << "I";
			 if (res == "LEU") output << "L";
			 if (res == "ARG") output << "R";
			 if (res == "LYS") output << "K";
			 if (res == "HIS") output << "H";
			 if (res == "TYR") output << "Y";
			 if (res == "TRP") output << "W";
			 if (res == "PHE") output << "F";
			 output << atoms[(allassis[best][j].atoms[0])].resnr << " ";
		}
		output << endl;
		scores[best] = -1;
	}
	output.close();
	output.clear();
	
	return;
}


void shiftassis(vector<vector<noe> > & allassis, vector<atom> atoms, vector<noe> noes, string csfile, string logfile){
	vector<string> data;
	vector<atom> shifts;
	vector<float> scores;
	ofstream output;
	
	data = readfile(csfile);
	for (unsigned int i = 0; i < data.size(); i++){
		stringstream ss;
		unsigned int i1;
		unsigned int i2;
		string s1;
		string s2;
		float f1;
		atom peak;
		
		ss << data[i];
		ss >> i1 >> i2 >> s1 >> s2 >> f1;
		if (s2 != "H") continue;
		peak.resnr = i2;
		peak.res = s1;
		peak.h = f1;
		shifts.push_back(peak);
	}


	for (unsigned int i = 0; i < data.size(); i++){
		stringstream ss;
		unsigned int i1;
		unsigned int i2;
		string s1;
		string s2;
		float f1;
		
		ss << data[i];
		ss >> i1 >> i2 >> s1 >> s2 >> f1;
		if (s2 != "N") continue;
		for (unsigned int j = 0; j < shifts.size(); j++) if (shifts[j].res == s1 && shifts[j].resnr == i2) shifts[j].n = f1;		
	}		
	
	
	//for (unsigned int i = 0; i < shifts.size(); i++){
	//	if (shifts[i].n < 108) shifts[i].n = shifts[i].n + 20.9;
	//	if (shifts[i].n > 130) shifts[i].n = shifts[i].n - 20.9;
	//}
	
	
	for (unsigned int i = 0; i < allassis.size(); i++){
		float score = 0;
		for (unsigned int j = 0; j < allassis[i].size(); j++) if (abs(allassis[i][j].h1 - allassis[i][j].h2) < 0.01 && abs(allassis[i][j].n1 - allassis[i][j].n2) < 0.1){
			float value = 0;
			for (unsigned int k = 0; k < shifts.size(); k++) if (shifts[k].resnr == atoms[(allassis[i][j].atoms[0])].resnr){
				// cout << allassis[i][j].n1 << " - " << shifts[k].n << endl;
				value = abs(allassis[i][j].n1 - shifts[k].n);
				value = value * 10;
				value = value + abs(allassis[i][j].h1 - shifts[k].h);
			}
			score = score + value;
		}
		scores.push_back(score);
	}
	
	for (unsigned int i = 0; i < scores.size(); i++) cout << fixed << setprecision(2) << scores[i] << endl;

	output.open(logfile.c_str(),ios::app);
	output << "Overview over ranked remaining assignments. Scored by the difference in the chemical shifts (hydrogen shifts times 10)" << endl;
	for (unsigned int i = 0; i < scores.size(); i++){	
		float score = -1;
		unsigned int best = 0;
		for (unsigned int j = 0; j < scores.size(); j++) if (scores[j] > 0) if (score < 0 || scores[j] < score){
			score = scores[j];
			best = j;
		}
		output << fixed << setprecision(1) << score << " ";
		for (unsigned int j = 0; j < allassis[best].size(); j++) if (abs(allassis[best][j].h1 - allassis[best][j].h2) < 0.01 && abs(allassis[best][j].n1 - allassis[best][j].n2) < 0.1){
			 string res = atoms[(allassis[best][j].atoms[0])].res;
			 if (res == "GLY") output << "G";
			 if (res == "ALA") output << "A";
			 if (res == "SER") output << "S";
			 if (res == "THR") output << "T";
			 if (res == "GLN") output << "Q";
			 if (res == "GLU") output << "E";
			 if (res == "ASN") output << "N";
			 if (res == "ASP") output << "D";
			 if (res == "CYS") output << "C";
			 if (res == "MET") output << "M";
			 if (res == "VAL") output << "V";
			 if (res == "ILE") output << "I";
			 if (res == "LEU") output << "L";
			 if (res == "ARG") output << "R";
			 if (res == "LYS") output << "K";
			 if (res == "HIS") output << "H";
			 if (res == "TYR") output << "Y";
			 if (res == "TRP") output << "W";
			 if (res == "PHE") output << "F";
			 output << atoms[(allassis[best][j].atoms[0])].resnr << " ";
		}
		output << endl;
		scores[best] = -1;
	}
	output.close();
	output.clear();


	
	return;
}

