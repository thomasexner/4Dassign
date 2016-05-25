unsigned int optimize(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, vector<atom> fix){
	bool ende = true;
	ofstream output;
	unsigned int errorlevel = 0;
	

	// Remove any random entries
	// Very important if called by artefacts routine!!
	for (unsigned int i = 0; i < noes.size(); i++){
		noes[i].atoms.clear();
		noes[i].atoms2.clear();
		noes[i].crosspeaks.clear();
		noes[i].assigned = false;
	}
   for (unsigned int i = 0; i < atoms.size(); i++){
		atoms[i].noes.clear();
		atoms[i].assigned = false;		
	}
	
	// Build crosspeaks vector	
	for (unsigned int i = 0; i < noes.size(); i++) for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1) noes[i].crosspeaks.push_back(j);
	// Writes all noes.atoms vectors: Create an inital assignment for every noe
	// Doublepeaks are not examined, they are just stored to the remaining spaces
	assignall(atoms, noes, safety, logfile, fix);
	
	// printassignment(atoms, noes, logfile);	
	printnoes(atoms, noes, logfile);
	
	// Check for holes: If there are holes return errorvalue(1) to main
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() == 0){
		cout << "Hole found in NOE: " << fixed << setprecision(2) << noes[i].h1 << " " << noes[i].n1 << " => " << noes[i].h2 << " " << noes[i].n2 << endl;
		return 1;
	}

	// Check if finished: If finished return 0
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() != 1) ende = false;
	if (ende == true) return 0;	
	
	// For the remaining not assigned NOEs go through all possibiblities, if one is violated remove it
	// Remove some entries and than recall reduceassignment routine
	cout << "Start optimizing routine" << endl;
	errorlevel = examineassignment(atoms, noes, logfile, safety);

	// printassignment(atoms, noes, logfile);
	printnoes(atoms, noes, logfile);
		
	return errorlevel;
}

unsigned int examineassignment(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety){
	// Go through the assignment possibilities and try each of them, the one that creates holes will be erased
	// Start with those noes that have the lowest number of entries
	while (1){
		unsigned int counts = 0;
		unsigned int newcounts = 0;
		bool ende = false;
		
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) counts = counts + noes[i].atoms.size();

		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].assigned == false){
			// cout << "Examine NOE: " << fixed << setprecision(2) << noes[i].n1 << " " << noes[i].h1 << " (" << noes[i].crosspeaks.size() << ") " << endl;
			if (noes[i].doublepeak == true){
				// Create a twodimensional matrix with the entries of the assignments
				// Use this matrix to check if an assignment is allow or not
				// If a comination of NOEs is not combinable with the current assignment vector label it as false, i.e. violated, in this matrix
				// Matrix will be used symmetric do not check assignments twice
				vector<vector<bool> > assmatrix;
				for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
					vector<bool> bools;
					for (unsigned int k = 0; k < noes[i].atoms.size(); k++) bools.push_back(true);
					assmatrix.push_back(bools);
				}
				for (unsigned int j = 0; j < assmatrix.size(); j++) assmatrix[j][j] = false;
			
				for (unsigned int j = 0; j < noes[i].atoms.size(); j++) for (unsigned int k = 0; k < noes[i].atoms.size(); k++) if (j < k){
					vector<noe> newnoes;
					vector<atom> newatoms;
					bool accept = true;
				
					for (unsigned int l = 0; l < noes.size(); l++) newnoes.push_back(noes[l]);
					for (unsigned int l = 0; l < atoms.size(); l++) newatoms.push_back(atoms[l]);
				
					newnoes[i].atoms.clear();
					newnoes[i].atoms.push_back(noes[i].atoms[j]);
					newnoes[i].atoms.push_back(noes[i].atoms[k]);
				
					reduceassignment(newatoms, newnoes, safety);
				
					accept = checkreduced(newatoms, newnoes, safety);
					
					if (accept == false){
						assmatrix[j][k] = false;
						assmatrix[k][j] = false;
					}
				}
				// Store all possibilities that are still allowed in keeps vector
				vector<unsigned int> keeps;
				for (unsigned int j = 0; j < assmatrix.size(); j++){
					bool accept = false;
					for (unsigned int k = 0; k < assmatrix[j].size();  k++) if (assmatrix[j][k] == true) accept = true;
					if (accept == true) keeps.push_back(noes[i].atoms[j]);
				}
				
				noes[i].atoms.clear();
				for (unsigned int j = 0; j < keeps.size(); j++){
					noes[i].atoms.push_back(keeps[j]);
				}
				reduceassignment(atoms, noes, safety);
			}
		
			if (noes[i].doublepeak == false) for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
				vector<noe> newnoes;
				vector<atom> newatoms;
				bool check = false;
				
				for (unsigned int k = 0; k < noes.size(); k++) newnoes.push_back(noes[k]);
				for (unsigned int k = 0; k < atoms.size(); k++) newatoms.push_back(atoms[k]);
		
				newnoes[i].atoms.clear();
				newnoes[i].atoms.push_back(noes[i].atoms[j]);
				
				reduceassignment(newatoms, newnoes, safety);
								
				check = checkreduced(newatoms, newnoes, safety);
				if (check == false){
					noes[i].atoms.erase(noes[i].atoms.begin() + j);
					j--;
				}
			}
		}
		reduceassignment(atoms, noes, safety);


		// Check if finished: If finished return 0
		ende = true;
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() != 1) ende = false;
		if (ende == true){
			cout << "Found final assignment. Return" << endl;
			return 0;	
		}

		
		// Check for holes, if holes give up assignment
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() == 0){
			// cout << "Hole found in NOE: " << fixed << setprecision(2) << noes[i].h1 << " " << noes[i].n1 << endl;
			return 1;
		}
				
		// printassignment(atoms, noes, logfile);
		// printnoes(atoms, noes, logfile);
		
		
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) newcounts = newcounts + noes[i].atoms.size();
		
		if (newcounts >= counts){
			// cout << "No assignment possibilities removed: " << counts << " => " << newcounts << " break" << endl;
			break;
		}

	}
	reduceassignment(atoms, noes, safety);

	return 0;
}

bool checkreduced(vector<atom> atoms, vector<noe> noes, float safety){
	// First check for holes, if there is a hole return false
	// otherwise go to second check routine below
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1){
		if (noes[i].atoms.size() == 0) return false;
		if (noes[i].atoms.size() == 1 && noes[i].doublepeak == true) return false; 
	}
	
	// Now check if any assigned noes conflicts with another asssigned noe
	// If any crosspeak between two assigned noes is violated assignment can be considered as violated
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n2 - noes[i].n1) < 0.1) if (noes[i].assigned == true) for (unsigned int j = 0; j < noes[i].crosspeaks.size(); j++){
		unsigned int diagonal = 0;
		for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h2 - noes[k].h1) < 0.01 && abs(noes[k].n1 - noes[k].n2) < 0.1){
			if (abs(noes[k].h1 - noes[(noes[i].crosspeaks[j])].h2) < 0.01 && abs(noes[k].n1 - noes[(noes[i].crosspeaks[j])].n2) < 0.1) diagonal = k;
		}
		if (noes[diagonal].assigned == false) continue;
		float dist = 0;
		float resi = 0;
		
		resi = noes[(noes[i].crosspeaks[j])].dist;
		resi = resi + safety;
		resi = resi * resi;
		
		dist = dist + ((atoms[(noes[i].atoms[0])].x - atoms[(noes[diagonal].atoms[0])].x) * (atoms[(noes[i].atoms[0])].x - atoms[(noes[diagonal].atoms[0])].x));
		dist = dist + ((atoms[(noes[i].atoms[0])].y - atoms[(noes[diagonal].atoms[0])].y) * (atoms[(noes[i].atoms[0])].y - atoms[(noes[diagonal].atoms[0])].y));
		dist = dist + ((atoms[(noes[i].atoms[0])].z - atoms[(noes[diagonal].atoms[0])].z) * (atoms[(noes[i].atoms[0])].z - atoms[(noes[diagonal].atoms[0])].z));
		
		if (noes[i].doublepeak == true){
			float newdist = 0;
			newdist = newdist + ((atoms[(noes[i].atoms[1])].x - atoms[(noes[diagonal].atoms[0])].x) * (atoms[(noes[i].atoms[1])].x - atoms[(noes[diagonal].atoms[0])].x));
			newdist = newdist + ((atoms[(noes[i].atoms[1])].y - atoms[(noes[diagonal].atoms[0])].y) * (atoms[(noes[i].atoms[1])].y - atoms[(noes[diagonal].atoms[0])].y));
			newdist = newdist + ((atoms[(noes[i].atoms[1])].z - atoms[(noes[diagonal].atoms[0])].z) * (atoms[(noes[i].atoms[1])].z - atoms[(noes[diagonal].atoms[0])].z));
			if (newdist < dist) dist = newdist;
		}

		if (noes[diagonal].doublepeak == true){
			float newdist = 0;
			newdist = newdist + ((atoms[(noes[i].atoms[0])].x - atoms[(noes[diagonal].atoms[1])].x) * (atoms[(noes[i].atoms[0])].x - atoms[(noes[diagonal].atoms[1])].x));
			newdist = newdist + ((atoms[(noes[i].atoms[0])].y - atoms[(noes[diagonal].atoms[1])].y) * (atoms[(noes[i].atoms[0])].y - atoms[(noes[diagonal].atoms[1])].y));
			newdist = newdist + ((atoms[(noes[i].atoms[0])].z - atoms[(noes[diagonal].atoms[1])].z) * (atoms[(noes[i].atoms[0])].z - atoms[(noes[diagonal].atoms[1])].z));
			if (newdist < dist) dist = newdist;
		}
		if (noes[i].doublepeak == true && noes[diagonal].doublepeak == true){
			float newdist = 0;
			newdist = newdist + ((atoms[(noes[i].atoms[1])].x - atoms[(noes[diagonal].atoms[1])].x) * (atoms[(noes[i].atoms[1])].x - atoms[(noes[diagonal].atoms[1])].x));
			newdist = newdist + ((atoms[(noes[i].atoms[1])].y - atoms[(noes[diagonal].atoms[1])].y) * (atoms[(noes[i].atoms[1])].y - atoms[(noes[diagonal].atoms[1])].y));
			newdist = newdist + ((atoms[(noes[i].atoms[1])].z - atoms[(noes[diagonal].atoms[1])].z) * (atoms[(noes[i].atoms[1])].z - atoms[(noes[diagonal].atoms[1])].z));
			if (newdist < dist) dist = newdist;
			
		}

		if (dist > resi) return false;
	}
	
	return true;
}


void reduceassignment(vector<atom> & atoms, vector<noe> & noes, float safety){
	// First check for assigned NOEs
	for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 1){
		noes[i].assigned = true;
		atoms[(noes[i].atoms[0])].assigned = true;
		atoms[(noes[i].atoms[0])].h = noes[i].h1;
		atoms[(noes[i].atoms[0])].n = noes[i].n1;
	}
	// Additionaly check for assigned NOE doublepeak
	for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 2 && noes[i].doublepeak == true){
		noes[i].assigned = true;
		atoms[(noes[i].atoms[0])].assigned = true;
		atoms[(noes[i].atoms[1])].assigned = true;
		atoms[(noes[i].atoms[0])].h = noes[i].h1;
		atoms[(noes[i].atoms[0])].n = noes[i].n1;
		atoms[(noes[i].atoms[1])].h = noes[i].h1;
		atoms[(noes[i].atoms[1])].n = noes[i].n1;
	}

	// Check for already assigned atoms
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n2 - noes[i].n1) < 0.1) if (noes[i].assigned == false) for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
		if (atoms[(noes[i].atoms[j])].assigned == true){
			noes[i].atoms.erase(noes[i].atoms.begin() + j);
			j--;
		}
	}

	// Reduce the assignments: Eliminate impossiblites due to created assignments
	while(1){
		unsigned int entries = 0;
		unsigned int newentries = 0;
		for (unsigned int i = 0; i < noes.size(); i++) entries = entries + noes[i].atoms.size();
		
		// Check if some entries can be eliminated
		// If any NOE has an atom that is assigned check if distance restraint is violated
		// If distance is violated remove this atom
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n2 - noes[i].n1) < 0.1) if (noes[i].doublepeak == false && noes[i].assigned == false) for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
			bool accept = true;
			
			for (unsigned int k = 0; k < noes[i].crosspeaks.size(); k++){
				// Now find corresponding diagonal peak
				for (unsigned int l = 0; l < noes.size(); l++) if (abs(noes[l].h1 - noes[l].h2) < 0.01 && abs(noes[l].n1 - noes[l].n2) < 0.1) if (abs(noes[l].h1 - noes[(noes[i].crosspeaks[k])].h2) < 0.01 && abs(noes[l].n1 - noes[(noes[i].crosspeaks[k])].n2) < 0.1){
					if (noes[l].assigned == true && noes[l].atoms.size() == 1){
						float dist = 0;
						float resi = noes[(noes[i].crosspeaks[k])].dist;
						
						resi = resi + safety;
						resi = resi * resi;
						
						dist = dist + ((atoms[(noes[l].atoms[0])].x - atoms[(noes[i].atoms[j])].x) * (atoms[(noes[l].atoms[0])].x - atoms[(noes[i].atoms[j])].x));
						dist = dist + ((atoms[(noes[l].atoms[0])].y - atoms[(noes[i].atoms[j])].y) * (atoms[(noes[l].atoms[0])].y - atoms[(noes[i].atoms[j])].y));
						dist = dist + ((atoms[(noes[l].atoms[0])].z - atoms[(noes[i].atoms[j])].z) * (atoms[(noes[l].atoms[0])].z - atoms[(noes[i].atoms[j])].z));
						
						if (dist > resi) accept = false;				
					}
					if (accept == false) break;
				}		
			}
			if (accept == false){
				noes[i].atoms.erase(noes[i].atoms.begin() + j);
				j--;
			}
		}
		
		// Check if there are any unqiue ones, if so remove this atom from other noe.atoms vectors
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() == 1){
			for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[j].h1 - noes[j].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) for (unsigned int k = 0; k < noes[j].atoms.size(); k++) if (noes[j].atoms[k] == noes[i].atoms[0]){
				noes[j].atoms.erase(noes[j].atoms.begin() + k);
				k--;
			}
		}
		
		// Check if any two noes share two atoms exclusively, if so remove both assignments from all other noes
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() == 2 && noes[i].doublepeak == false){
			for (unsigned int j = 0; j < noes.size(); j++) if (j != i && noes[j].doublepeak == false) if (abs(noes[j].h1 - noes[j].h2) < 0.01 && abs(noes[j].n2 - noes[j].n1) < 0.1) if (noes[j].atoms.size() == 2){
				bool found = false;
				if (noes[i].atoms[0] == noes[j].atoms[0] && noes[i].atoms[1] == noes[j].atoms[1]) found = true;
				if (noes[i].atoms[0] == noes[j].atoms[1] && noes[i].atoms[1] == noes[j].atoms[0]) found = true;
				if (found == false) continue;
				for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h1 - noes[k].h2) < 0.01 && abs(noes[k].n1 - noes[k].n2) < 0.1) if (k != i) if (k != j) for (unsigned int l = 0; l < noes[k].atoms.size(); l++){
					if (noes[k].atoms[l] == noes[i].atoms[0] || noes[k].atoms[l] == noes[i].atoms[1]){
						noes[k].atoms.erase(noes[k].atoms.begin() + l);
						l--;
					}
				}
			}
		}
		
		for (unsigned int i = 0; i < noes.size(); i++) newentries = newentries + noes[i].atoms.size();
		if (newentries < entries){
			// Label assigned peaks and store shifts in atoms vector
			for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 1){
				noes[i].assigned = true;
				atoms[(noes[i].atoms[0])].assigned = true;
				atoms[(noes[i].atoms[0])].h = noes[i].h1;
				atoms[(noes[i].atoms[0])].n = noes[i].n1;
			}

			for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 2 && noes[i].doublepeak == true){
				noes[i].assigned = true;
				atoms[(noes[i].atoms[0])].assigned = true;
				atoms[(noes[i].atoms[1])].assigned = true;
				atoms[(noes[i].atoms[0])].h = noes[i].h1;
				atoms[(noes[i].atoms[0])].n = noes[i].n1;
				atoms[(noes[i].atoms[1])].h = noes[i].h1;
				atoms[(noes[i].atoms[1])].n = noes[i].n1;
			}
			continue;
		}
		
		/*
		// Check triangular inequality: If B is assigned the distance A to C must be smaller or equal than the sum of the distances A to B plus B to C
		// If violated the assignment is impossible, so create a hole and return
		for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].assigned == true && noes[i].crosspeaks.size() > 1) if (abs(noes[i].n1 - noes[i].n2) < 0.1 && abs(noes[i].h1 - noes[i].h2) < 0.01) for (unsigned int j = 0; j < noes[i].crosspeaks.size(); j++){
			float dist1 = noes[(noes[i].crosspeaks[j])].dist;
			for (unsigned int k = 0; k < noes[i].crosspeaks.size(); k++){
				if (k == j) continue;
				float dist2 = noes[(noes[i].crosspeaks[k])].dist;
				unsigned int peakpos = 0;
				bool found = false;
				
				for (unsigned int m = 0; m < noes.size(); m++){
					if (abs(noes[m].h1 - noes[(noes[i].crosspeaks[j])].h2) < 0.01 && abs(noes[m].n1 - noes[(noes[i].crosspeaks[j])].n2) < 0.1) if (abs(noes[m].h2 - noes[(noes[i].crosspeaks[k])].h2) < 0.01 && abs(noes[m].n2 - noes[(noes[i].crosspeaks[k])].n2) < 0.1){
						found = true;
						peakpos = m;
						break;
					}
					if (abs(noes[m].h2 - noes[(noes[i].crosspeaks[j])].h2) < 0.01 && abs(noes[m].n2 - noes[(noes[i].crosspeaks[j])].n2) < 0.1) if (abs(noes[m].h1 - noes[(noes[i].crosspeaks[k])].h2) < 0.01 && abs(noes[m].n1 - noes[(noes[i].crosspeaks[k])].n2) < 0.1){
						found = true;
						peakpos = m;
						break;
					}
				}
				if (found == false) continue;

				// Found a triangular equation, if violated create a hole in i and return to previous routine
				if (noes[peakpos].dist > (dist1 + dist2)){
					cout << "Remove TE: " << noes[(noes[i].crosspeaks[j])].h1 << " " << noes[(noes[i].crosspeaks[j])].n1 << " => " << noes[(noes[i].crosspeaks[j])].h2 << " " << noes[(noes[i].crosspeaks[j])].n2 << endl;
					noes[i].assigned = false;
					noes[i].atoms.clear();
					return;
				
				}
				
			}
			
		}
		*/
		
		
		// Check if for any assignment entry all combinations are violated
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].assigned == false && noes[i].doublepeak == false) for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
			bool accept = true;
			for (unsigned int k = 0; k < noes[i].crosspeaks.size(); k++){
				// Find the corresponding diagonal peak
				// Check if the restraint is violated for all assignments of this peak
				for (unsigned int l = 0; l < noes.size(); l++) if (abs(noes[l].h1 - noes[l].h2) < 0.01 && abs(noes[l].n2 - noes[l].n1) < 0.1) if (abs(noes[l].h1 - noes[(noes[i].crosspeaks[k])].h2) < 0.01 && abs(noes[l].n1 - noes[(noes[i].crosspeaks[k])].n2) < 0.1){
					if (noes[l].assigned == true) break;
					bool check = false;
					// If there is at least one possibility to fullfill the restraint keep the assignment (label accept still true)
					for (unsigned int m = 0; m < noes[l].atoms.size(); m++){
						float dist = 0;
						float resi = noes[(noes[i].crosspeaks[k])].dist;
						
						resi = resi + safety;
						resi = resi * resi;
						
						dist = dist + ((atoms[(noes[i].atoms[j])].x - atoms[(noes[l].atoms[m])].x) * (atoms[(noes[i].atoms[j])].x - atoms[(noes[l].atoms[m])].x));
						dist = dist + ((atoms[(noes[i].atoms[j])].y - atoms[(noes[l].atoms[m])].y) * (atoms[(noes[i].atoms[j])].y - atoms[(noes[l].atoms[m])].y));
						dist = dist + ((atoms[(noes[i].atoms[j])].z - atoms[(noes[l].atoms[m])].z) * (atoms[(noes[i].atoms[j])].z - atoms[(noes[l].atoms[m])].z));
						
						if (dist < resi) check = true;
					}
					if (check == false) accept = false;
				}
				// If there is already a violation do not check the rest but remove
				if (accept == false) break;
			}
			
			if (accept == false){
				noes[i].atoms.erase(noes[i].atoms.begin() + j);
				j--;
			}
		}
		
		for (unsigned int i = 0; i < noes.size(); i++) newentries = newentries + noes[i].atoms.size();
		if (newentries < entries){
			// Label assigned peaks and store shifts in atoms vector
			for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 1){
				noes[i].assigned = true;
				atoms[(noes[i].atoms[0])].assigned = true;
				atoms[(noes[i].atoms[0])].h = noes[i].h1;
				atoms[(noes[i].atoms[0])].n = noes[i].n1;
			}

			for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 2 && noes[i].doublepeak == true){
				noes[i].assigned = true;
				atoms[(noes[i].atoms[0])].assigned = true;
				atoms[(noes[i].atoms[1])].assigned = true;
				atoms[(noes[i].atoms[0])].h = noes[i].h1;
				atoms[(noes[i].atoms[0])].n = noes[i].n1;
				atoms[(noes[i].atoms[1])].h = noes[i].h1;
				atoms[(noes[i].atoms[1])].n = noes[i].n1;
			}
			continue;
		}
	
		
		// Now the doublepeak routines:
		// Now check if any crosspeak points on an already assigned atom, if so check if distance restraint is violated,
		// if so remove the atoms vector, there is no assignment possible for this vector
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n2 - noes[i].n1) < 0.1) if (noes[i].doublepeak == true && noes[i].assigned == false) for (unsigned int j = 0; j < noes[i].crosspeaks.size(); j++){
			// Find the corresponding diagonalpeak
			unsigned int diagonal = 0;
			float dist = 0;
			float resi = 0;
			for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h1 - noes[k].h2) < 0.01 && abs(noes[k].n1 - noes[k].n1) < 0.1){
				if (abs(noes[k].h1 - noes[(noes[i].crosspeaks[j])].h2) < 0.01 && abs(noes[k].n1 - noes[(noes[i].crosspeaks[j])].n2) < 0.1) diagonal = k;
			}
			if (noes[diagonal].assigned == false) continue;
			for (unsigned int k = 0; k < noes[i].atoms.size(); k++){
				float newdist = 0;		
				
				newdist = newdist + ((atoms[(noes[i].atoms[k])].x - atoms[(noes[diagonal].atoms[0])].x) * (atoms[(noes[i].atoms[k])].x - atoms[(noes[diagonal].atoms[0])].x));
				newdist = newdist + ((atoms[(noes[i].atoms[k])].y - atoms[(noes[diagonal].atoms[0])].y) * (atoms[(noes[i].atoms[k])].y - atoms[(noes[diagonal].atoms[0])].y));
				newdist = newdist + ((atoms[(noes[i].atoms[k])].z - atoms[(noes[diagonal].atoms[0])].z) * (atoms[(noes[i].atoms[k])].z - atoms[(noes[diagonal].atoms[0])].z));
				if (dist < 0 || newdist < dist) dist = newdist;
			}
						
			resi = noes[(noes[i].crosspeaks[j])].dist;
			resi = resi + safety;
			resi = resi * resi;

			if (dist > resi) noes[i].atoms.clear();
		}	 
		// Now check if there is an assigned doublepeak 
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].doublepeak == true && noes[i].assigned == true) for (unsigned int j = 0; j < noes[i].crosspeaks.size(); j++){
			if (noes[i].atoms.size() != 2) continue;
			unsigned int diagonal = 0;
			for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h1 - noes[k].h2) < 0.01 && abs(noes[k].n1 - noes[k].n2) < 0.1){
				if (abs(noes[k].h1 - noes[(noes[i].crosspeaks[j])].h2) < 0.01 && abs(noes[k].n1 - noes[(noes[i].crosspeaks[j])].n2) < 0.1) diagonal = k;
			}
			if (noes[diagonal].assigned == true) continue;
			for (unsigned int k = 0; k < noes[diagonal].atoms.size(); k++){
				float dist = 0;
				float newdist = 0;
				float resi = 0;
				
				dist = dist + ((atoms[(noes[diagonal].atoms[k])].x - atoms[(noes[i].atoms[0])].x) * (atoms[(noes[diagonal].atoms[k])].x - atoms[(noes[i].atoms[0])].x));
				dist = dist + ((atoms[(noes[diagonal].atoms[k])].y - atoms[(noes[i].atoms[0])].y) * (atoms[(noes[diagonal].atoms[k])].y - atoms[(noes[i].atoms[0])].y));
				dist = dist + ((atoms[(noes[diagonal].atoms[k])].z - atoms[(noes[i].atoms[0])].z) * (atoms[(noes[diagonal].atoms[k])].z - atoms[(noes[i].atoms[0])].z));

				newdist = newdist + ((atoms[(noes[diagonal].atoms[k])].x - atoms[(noes[i].atoms[1])].x) * (atoms[(noes[diagonal].atoms[k])].x - atoms[(noes[i].atoms[1])].x));
				newdist = newdist + ((atoms[(noes[diagonal].atoms[k])].y - atoms[(noes[i].atoms[1])].y) * (atoms[(noes[diagonal].atoms[k])].y - atoms[(noes[i].atoms[1])].y));
				newdist = newdist + ((atoms[(noes[diagonal].atoms[k])].z - atoms[(noes[i].atoms[1])].z) * (atoms[(noes[diagonal].atoms[k])].z - atoms[(noes[i].atoms[1])].z));
				
				if (newdist < dist) dist = newdist;
				resi = noes[(noes[i].crosspeaks[j])].dist;
				resi = resi + safety;
				resi = resi * resi;
				if (dist > resi){
					noes[diagonal].atoms.erase(noes[diagonal].atoms.begin() + k);
					k--;
				}
			
			}
		}
				
		// Label assigned peaks and store shifts in atoms vector
		for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 1){
			noes[i].assigned = true;
			atoms[(noes[i].atoms[0])].assigned = true;
			atoms[(noes[i].atoms[0])].h = noes[i].h1;
			atoms[(noes[i].atoms[0])].n = noes[i].n1;
		}

		for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 2 && noes[i].doublepeak == true){
			noes[i].assigned = true;
			atoms[(noes[i].atoms[0])].assigned = true;
			atoms[(noes[i].atoms[1])].assigned = true;
			atoms[(noes[i].atoms[0])].h = noes[i].h1;
			atoms[(noes[i].atoms[0])].n = noes[i].n1;
			atoms[(noes[i].atoms[1])].h = noes[i].h1;
			atoms[(noes[i].atoms[1])].n = noes[i].n1;
		}

		// Check for already assigned atoms
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n2 - noes[i].n1) < 0.1) if (noes[i].assigned == false) for (unsigned int j = 0; j < noes[i].atoms.size(); j++){
			if (atoms[(noes[i].atoms[j])].assigned == true){
				noes[i].atoms.erase(noes[i].atoms.begin() + j);
				j--;
			}
		
		}
		
		
		for (unsigned int i = 0; i < noes.size(); i++) newentries = newentries + noes[i].atoms.size();
		if (newentries >= entries) break;
		
		// Check for holes:
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h2 - noes[i].h1) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1){
			if (noes[i].atoms.size() == 0) return;
			if (noes[i].doublepeak == true && noes[i].atoms.size() == 1) return;
		}
	}
	
	return;
}


void fix_atoms(vector<atom> & atoms, vector<noe> & noes, vector<atom> fix){
	cout << "Fix " << fix.size() << " NOESY-signals" << endl;
	if (fix.size() == 0) return;
	for (unsigned int i = 0; i < fix.size(); i++) for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[j].h1 - noes[j].h2) < 0.01 && abs(noes[j].n1 - noes[j].n2) < 0.1) if (abs(noes[j].h1 - fix[i].h) < 0.01 && abs(noes[j].n1 - fix[i].n) < 0.1){
		for (unsigned int k = 0; k < atoms.size(); k++) if (atoms[k].resnr == fix[i].resnr){
			noes[j].atoms.clear();
			noes[j].atoms.push_back(k);
			noes[i].assigned = true;
			atoms[k].assigned = true;
			atoms[k].h = noes[j].h1;
			atoms[k].n = noes[j].n1;
		}
	}
	
	return;
}

void assignall(vector<atom> & atoms, vector<noe> & noes, float safety, string logfile, vector<atom> fix){
	// Create an initial assignment
	cout << "Get initial assignment for " << noes.size() << " peaks." << endl;
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1){
		if (noes[i].nh2 == true) continue;
		if (noes[i].doublepeak == true) continue;
		
		for (unsigned int j = 0; j < atoms.size(); j++) if (atoms[j].assigned == false){
			if (atoms[j].name == "HE21" && atoms[j].res == "GLN") continue;
			if (atoms[j].name == "HE22" && atoms[j].res == "GLN") continue;
			if (atoms[j].name == "HD21" && atoms[j].res == "ASN") continue;
			if (atoms[j].name == "HD22" && atoms[j].res == "ASN") continue;
			bool check = true;
			if (noes[i].crosspeaks.size() < 7) check = checkassignment(atoms, noes, j, i, safety);
			if (check == true){
				noes[i].atoms.push_back(j);	
				atoms[j].noes.push_back(i);
			}
			
		}
	}
	// cout << "Call speedup routine" << endl;

	// printnoes(atoms, noes, logfile);
		
	fix_atoms(atoms, noes, fix);
	
	// speedup(noes, atoms);
	
	// printnoes(atoms, noes, logfile);
		
	// Reduce assignment for the first time get some assigned peaks:
	while(1){
		unsigned int assis = 0;
		unsigned int newassis = 0;
		for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].assigned == true) assis++;
		
		// Check if there are any unqiue ones, if so remove this atom from other noe.atoms vectors
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() == 1){
			for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[j].h1 - noes[j].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) for (unsigned int k = 0; k < noes[j].atoms.size(); k++) if (noes[j].atoms[k] == noes[i].atoms[0]){
				noes[j].atoms.erase(noes[j].atoms.begin() + k);
				k--;
			}
			
		}
		// Check if any two noes share two atoms exclusively, if so remove both assignments from all other noes
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].atoms.size() == 2){
			for (unsigned int j = 0; j < noes.size(); j++) if (j != i) if (abs(noes[j].h1 - noes[j].h2) < 0.01 && abs(noes[j].n2 - noes[j].n1) < 0.1) if (noes[j].atoms.size() == 2){
				bool found = false;
				if (noes[i].atoms[0] == noes[j].atoms[0] && noes[i].atoms[1] == noes[j].atoms[1]) found = true;
				if (noes[i].atoms[0] == noes[j].atoms[1] && noes[i].atoms[1] == noes[j].atoms[0]) found = true;
				if (found == false) continue;
				
				for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h1 - noes[k].h2) < 0.01 && abs(noes[k].n1 - noes[k].n2) < 0.1) if (k != i) if (k != j) for (unsigned int l = 0; l < noes[k].atoms.size(); l++){
					if (noes[k].atoms[l] == noes[i].atoms[0] || noes[k].atoms[l] == noes[i].atoms[1]){
						noes[k].atoms.erase(noes[k].atoms.begin() + l);
						l--;
					}
				}
			}
		}
		// Label assigned peaks and store shifts in atoms vector
		for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].atoms.size() == 1){
			noes[i].assigned = true;
			atoms[(noes[i].atoms[0])].assigned = true;
			atoms[(noes[i].atoms[0])].h = noes[i].h1;
			atoms[(noes[i].atoms[0])].n = noes[i].n1;
		}
	
		// If no new assignments are made break the loop
		for (unsigned int i = 0; i < noes.size(); i++) if (noes[i].assigned == true) newassis++;
		if (newassis <= assis) break;
	}
	// printnoes(atoms, noes, logfile);
	
	// Now add the doublepeak(s): Usually only one or a few peaks in the spectrum are doublepeaks
	// Assume that there will be no doublepeaks of NH2 sidechains
	// Full assignment check is to time consuming, so store all non assigned backbones here:
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].doublepeak == true) for (unsigned int j = 0; j < atoms.size(); j++){
		if (atoms[j].res == "ASN" && atoms[j].name == "HD21") continue;
		if (atoms[j].res == "ASN" && atoms[j].name == "HD22") continue;
		if (atoms[j].res == "GLN" && atoms[j].name == "HE21") continue;
		if (atoms[j].res == "GLN" && atoms[j].name == "HE22") continue;
		if (atoms[j].assigned == true) continue;
		noes[i].atoms.push_back(j);
		noes[i].atoms2.push_back(j);
		atoms[j].noes.push_back(i);
	}
	// printnoes(atoms, noes, logfile);
	// Reduce assignment: Now take the already assigned atoms into account and remove those entries from the noes.atoms vector that have assigned atoms that violate the distance restraint
	reduceassignment(atoms, noes, safety);
	// printnoes(atoms, noes, logfile);
	return;
}

bool checkassignment(vector<atom> atoms, vector<noe> noes, unsigned int atompeak, unsigned int noepeak, float safety){
	bool check = false;
	vector<vector<unsigned int> > assignment;
	vector<unsigned int> assis;
	
	// Create an assignment vector which has an entry for every crosspeak: 
	// Each of these entries is a vector which contains all valid assignments for this crosspeak
	// Find for each crosspeak every H(N) atom that fullfills the distance restraint 
	// Store these ones in the corresponding line of assignment
	for (unsigned int i = 0; i < noes[noepeak].crosspeaks.size(); i++){
		assis.clear();
		for (unsigned int j = 0; j < atoms.size(); j++){
			if (j == atompeak) continue;
			float dist = 0;
			float resi = noes[(noes[(noepeak)].crosspeaks[i])].dist;
			bool atomnh2 = false;
			bool noenh2 = false;
			
			// Find correspondig diagonalpeak on which this crosspeak is pointing
			// If this is an NH2 peak and the atom not (or vice versa) do not make an assignment
			for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h2 - noes[k].h1) < 0.01 && abs(noes[k].n2 - noes[k].n1) < 0.1){
				if (abs(noes[k].n1 - noes[(noes[noepeak].crosspeaks[i])].n2) < 0.1 && abs(noes[k].h1 - noes[(noes[noepeak].crosspeaks[i])].h2) < 0.01) if (noes[k].nh2 == true) noenh2 = true;
			}
			if (atoms[j].name == "HE21" && atoms[j].res == "GLN") atomnh2 = true;
			if (atoms[j].name == "HE22" && atoms[j].res == "GLN") atomnh2 = true;
			if (atoms[j].name == "HD21" && atoms[j].res == "ASN") atomnh2 = true;
			if (atoms[j].name == "HD22" && atoms[j].res == "ASN") atomnh2 = true;
			
			if (atomnh2 == true && noenh2 == false) continue;
			if (atomnh2 == false && noenh2 == true) continue;
			
			
			resi = resi + safety;
			resi = resi * resi;
			dist = dist + ((atoms[atompeak].x - atoms[j].x) * (atoms[atompeak].x - atoms[j].x));
			dist = dist + ((atoms[atompeak].y - atoms[j].y) * (atoms[atompeak].y - atoms[j].y));
			dist = dist + ((atoms[atompeak].z - atoms[j].z) * (atoms[atompeak].z - atoms[j].z));
			if (dist < resi) assis.push_back(j);
		}
		assignment.push_back(assis);
	}

	
	// First check if there are holes, if this is the case it is impossible to find an assignment for this crosspeak and the assignment check failed
	for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 0) return false;
	
	
	// Check if there are any unique ones: This means if for one crosspeak there is only one atom that fullfilles the distance restraint,
	// do not allow that this atom is assigned to any other crosspeak
	// Additionally check if there are any pairs, two peaks share to values, if so remove this from all others 
	while(1){
		unsigned int entries = 0;
		unsigned int newentries = 0;
		bool unique;
		for (unsigned int i = 0; i < assignment.size(); i++) entries = entries + assignment[i].size();

		// If there are only uniques end procedure
		unique = true;
		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() > 1) unique = false;
		if (unique == true) break; 
		
		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 1){
			for (unsigned int j = 0; j < assignment.size(); j++) if (j != i) for (unsigned int k = 0; k < assignment[j].size(); k++) if (assignment[j][k] == assignment[i][0]){
				assignment[j].erase(assignment[j].begin() + k);
				k--;
			}
		}

		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 2) for (unsigned int j = 0; j < assignment.size(); j++) if (j != i) if (assignment[j].size() == 2){
			bool found = false;
			if (assignment[i][0] == assignment[j][0] && assignment[j][1] == assignment[i][1]) found = true;
			if (assignment[i][0] == assignment[j][1] && assignment[i][1] == assignment[j][0]) found = true;
			if (found == false) continue;
			unsigned int i1 = assignment[i][0];
			unsigned int i2 = assignment[i][1];
			assignment[i].clear();
			assignment[j].clear();
			assignment[i].push_back(i1);
			assignment[j].push_back(i2);
		}		

		for (unsigned int i = 0; i < assignment.size(); i++) newentries = newentries + assignment[i].size();
		if (newentries >= entries) break;

		// If any hole is created break loop and return false
		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 0)	break;
		
	}

	// Check again for holes in the assignment
	for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 0)	return false;

	

	// Going recursively through all assignment possibilities: Pick a number from each line and store them in a vector:
	// If this vector has only different entries it is a valid assignment since now for every crosspeak a different atom can be assigned
	// If for at least one permutation such a valid assignment can be found label check as true
	// check is a pointer (always the same value) while assignment is used differently depending on the layer (subcall) of the routine
	check = false;
	getmatrix(assignment, check);
	
	
	return check;
}


bool checkout(vector<vector<unsigned int> > assignment){
	vector<unsigned int> assis;
	for (unsigned int i = 0; i < assignment.size(); i++) assis.push_back(assignment[i][0]);
	
	bool check = true;
	for (unsigned int i1 = 0; i1 < assis.size(); i1++) for (unsigned int i2 = 0; i2 < assis.size(); i2++) if (i2 != i1) if (assis[i2] == assis[i1]) check = false;

	return check;		
		
}


void getmatrix(vector<vector<unsigned int> > assignment, bool & check){
	// First check if the current assignment matrix is unique, i.e. each line has exactly one entry, if so call checking routine 'checkout' to verify if assignment is valid
	bool finished = true;
	for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() != 1) finished = false;
	
	
	
	if (finished == true){
		// Only overwrite check if it is true (check is initialized with false)
		if (checkout(assignment) == true) check = true;
		return;
	}
	
	// Go through assignment matrix: find the first line that has more than one entry
	// Create out of that line as many submatrices (newassignment) that have one entry in this line as there are entries in that line
	// Call this getmatrix routine with each of this submatrices
	// This will end up with (many) vectors that will call the checking routine above
	for (unsigned int i = 0; i < assignment.size() && check == false; i++){
		if (assignment[i].size() == 1) continue;
		for (unsigned int j = 0; j < assignment[i].size() && check == false; j++){
			bool accept = true;
			for (unsigned int k = 0; k < assignment.size(); k++) if (assignment[k].size() == 1 && k < j) if (assignment[k][0] == assignment[i][j]) accept = false;
			if (accept == false) continue;
			vector<vector<unsigned int> > newassignment;
			for (unsigned int k = 0; k < assignment.size(); k++) newassignment.push_back(assignment[k]);
			newassignment[i].clear();
			newassignment[i].push_back(assignment[i][j]);
			getmatrix(newassignment, check);
		}
		break;
	}
	
	
	return;
}

bool checkdoublepeak(vector<atom> atoms, vector<noe> noes, unsigned int noe0, unsigned int noe1, unsigned int noe2, unsigned int atom00, unsigned int atom01, unsigned int atom10, unsigned int atom11, float safety){
	bool check = false;
	
	// Only noe0 is a doublepeak: So first check if noe1 can be assigned to atom01 and noe2 can be assigned to atom11
	check = checkassignment(atoms, noes, atom01, noe1, safety);
	if (check == false) return check;
	
	check = checkassignment(atoms, noes, atom11, noe2, safety);
	if (check == false) return check;
	
	// If arrived here atom01 and atom11 can be assigned to noe1 and noe2 respectively so check what the doublepeak does:
	// Create one matrix for both atoms:
	// Use routine similar to above mentioned checkassignment one, except for creating the assignment matrix
	vector<vector<unsigned int> > assignment;
	vector<unsigned int> assis;

	// Create an assignment vector which has an entry for every crosspeak: 
	// Each of these entries is a vector which contains all valid assignments for this crosspeak
	// Find for each crosspeak every H(N) atom that fullfills the distance restraint 
	// Store these ones in the corresponding line of assignment
	// Before storing them, check if it is
	for (unsigned int i = 0; i < noes[noe0].crosspeaks.size(); i++){
		assis.clear();
		for (unsigned int j = 0; j < atoms.size(); j++){
			float dist = 0;
			float resi = noes[(noes[(noe0)].crosspeaks[i])].dist;
			bool atomnh2 = false;
			bool noenh2 = false;
			bool accept = false;
			
			// Find correspondig diagonalpeak on which this crosspeak is pointing
			// If this is an NH2 peak and the atom not (or vice versa) do not make an assignment
			for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[k].h2 - noes[k].h1) < 0.01 && abs(noes[k].n2 - noes[k].n1) < 0.1){
				if (abs(noes[k].n1 - noes[(noes[noe0].crosspeaks[i])].n2) < 0.1 && abs(noes[k].h1 - noes[(noes[noe0].crosspeaks[i])].h2) < 0.01) if (noes[k].nh2 == true) noenh2 = true;
			}
			if (atoms[j].name == "HE21" && atoms[j].res == "GLN") atomnh2 = true;
			if (atoms[j].name == "HE22" && atoms[j].res == "GLN") atomnh2 = true;
			if (atoms[j].name == "HD21" && atoms[j].res == "ASN") atomnh2 = true;
			if (atoms[j].name == "HD22" && atoms[j].res == "ASN") atomnh2 = true;
			
			if (atomnh2 == true && noenh2 == false) continue;
			if (atomnh2 == false && noenh2 == true) continue;
			
			
			resi = resi + safety;
			resi = resi * resi;
			if (j != atom00){
				dist = dist + ((atoms[atom00].x - atoms[j].x) * (atoms[atom00].x - atoms[j].x));
				dist = dist + ((atoms[atom00].y - atoms[j].y) * (atoms[atom00].y - atoms[j].y));
				dist = dist + ((atoms[atom00].z - atoms[j].z) * (atoms[atom00].z - atoms[j].z));
				if (dist < resi) accept = true;
			}
			if (j != atom10 && accept == false){
				dist = 0;
				dist = dist + ((atoms[atom10].x - atoms[j].x) * (atoms[atom10].x - atoms[j].x));
				dist = dist + ((atoms[atom10].y - atoms[j].y) * (atoms[atom10].y - atoms[j].y));
				dist = dist + ((atoms[atom10].z - atoms[j].z) * (atoms[atom10].z - atoms[j].z));
				if (dist < resi) accept = true;
			}
			if (accept == true) assis.push_back(j);
		}
		assignment.push_back(assis);
	}
	
	
	// First check if there are holes, if this is the case it is impossible to find an assignment for this crosspeak and the assignment check failed
	for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 0) return false;
	
	
	// Check if there are any unique ones: This means if for one crosspeak there is only one atom that fullfilles the distance restraint,
	// do not allow that this atom is assigned to any other crosspeak
	// Additionally check if there are any pairs, two peaks share to values, if so remove this from all others 
	while(1){
		unsigned int entries = 0;
		unsigned int newentries = 0;
		bool unique;
		for (unsigned int i = 0; i < assignment.size(); i++) entries = entries + assignment[i].size();

		// If there are only uniques end procedure
		unique = true;
		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() > 1) unique = false;
		if (unique == true) break; 
		
		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 1){
			for (unsigned int j = 0; j < assignment.size(); j++) if (j != i) for (unsigned int k = 0; k < assignment[j].size(); k++) if (assignment[j][k] == assignment[i][0]){
				assignment[j].erase(assignment[j].begin() + k);
				k--;
			}
		}

		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 2) for (unsigned int j = 0; j < assignment.size(); j++) if (j != i) if (assignment[j].size() == 2){
			bool found = false;
			if (assignment[i][0] == assignment[j][0] && assignment[j][1] == assignment[i][1]) found = true;
			if (assignment[i][0] == assignment[j][1] && assignment[i][1] == assignment[j][0]) found = true;
			if (found == false) continue;
			unsigned int i1 = assignment[i][0];
			unsigned int i2 = assignment[i][1];
			assignment[i].clear();
			assignment[j].clear();
			assignment[i].push_back(i1);
			assignment[j].push_back(i2);
		}		

		for (unsigned int i = 0; i < assignment.size(); i++) newentries = newentries + assignment[i].size();
		if (newentries >= entries) break;

		// If any hole is created break loop and return false
		for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 0)	break;
		
	}

	// Check again for holes in the assignment
	for (unsigned int i = 0; i < assignment.size(); i++) if (assignment[i].size() == 0)	return false;


	// Going recursively through all assignment possibilities: Pick a number from each line and store them in a vector:
	// If this vector has only different entries it is a valid assignment since now for every crosspeak a different atom can be assigned
	// If for at least one permutation such a valid assignment can be found label check as true
	// check is a pointer (always the same value) while assignment is used differently depending on the layer (subcall) of the routine
	check = false;
	getmatrix(assignment, check);
	
	
	return check;
}


unsigned int artefacts(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, vector<atom> fix){
	vector<noe> block;
	vector<noe> diags;
	vector<noe> cross;

	srand(time(NULL));
	cout << "Start artefact loop" << endl;
	while(1){	
		unsigned int result = 0;
		block.clear();
		diags.clear();
		cross.clear();
		for (unsigned int i = 0; i < noes.size(); i++){
			if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) diags.push_back(noes[i]);
			else cross.push_back(noes[i]);
		}
		for (unsigned int i = 0; i < diags.size(); i++) block.push_back(diags[i]);
		
		
		while(1){
			unsigned int guess = 0;
			unsigned int value = cross.size();
			
			// Guess a few random diagonal signals that have noe artefacts hopefully 
			while(cross.size() > block.size()){
				guess = rand() % cross.size();
				block.push_back(cross[guess]);
				cross.erase(cross.begin() + guess);
			}
			if (value > 100) for (unsigned int i = 0; i < 10; i++){
				guess = rand() % cross.size();
				block.push_back(cross[guess]);
				cross.erase(cross.begin() + guess);
			}
			else if (value > 50) for (unsigned int i = 0; i < 5; i++){
				guess = rand() % cross.size();
				block.push_back(cross[guess]);
				cross.erase(cross.begin() + guess);
			}
			else if (value > 20) for (unsigned int i = 0; i < 2; i++){
				guess = rand() % cross.size();
				block.push_back(cross[guess]);
				cross.erase(cross.begin() + guess);
			}
			else {
				guess = rand() % cross.size();
				block.push_back(cross[guess]);
				cross.erase(cross.begin() + guess);
			}
				
	
			cout << "Start block optimization with " << (block.size() - diags.size()) << " restraints" << endl;
			
			result = optimize(atoms,block,logfile,safety, fix);
			if (result != 0) break;
			
			if (result == 0){
				bool check = false;
				for (unsigned int i = 0; i < block.size(); i++) if (abs(block[i].h1 - block[i].h2) < 0.01 && abs(block[i].n1 - block[i].n2) < 0.1) if (block[i].assigned == true) check = true;
				if (check == true) printnoes(atoms, block, logfile);
			}
			
			// Check if there are any assigned peaks
			// If so check if there are any crosspeaks that are not included in the block are violated
			for (unsigned int i = 0; i < block.size(); i++) if (abs(block[i].h1 - block[i].h2) < 0.01 && abs(block[i].n1 - block[i].n2) < 0.1) if (block[i].assigned == true){
				for (unsigned int j = 0; j < block.size(); j++) if (j != i) if (abs(block[j].h1 - block[j].h2) < 0.01 && abs(block[j].n1 - block[j].n2) < 0.1) if (block[j].assigned == true){
					for (unsigned int k = 0; k < cross.size(); k++) if (abs(cross[k].h1 - block[i].h1) < 0.01 && abs(cross[k].n1 - block[i].n1) < 0.1) if (abs(cross[k].h2 - block[j].h2) < 0.01 && abs(cross[k].n2 - block[j].n2) < 0.1){
						float dist = 0;
						float restraint = 0;
						restraint = cross[k].dist + safety;
						restraint = restraint * restraint;
						
						dist = dist + ((atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[0])].x) * (atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[0])].x));  
						dist = dist + ((atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[0])].y) * (atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[0])].y));  
						dist = dist + ((atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[0])].z) * (atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[0])].z));  
						
						if (block[i].doublepeak == true){
							float newdist = 0;						
							newdist = newdist + ((atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[0])].x) * (atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[0])].x));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[0])].y) * (atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[0])].y));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[0])].z) * (atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[0])].z));  
							if (newdist < dist) dist = newdist;
						}
						if (block[j].doublepeak == true){
							float newdist = 0;						
							newdist = newdist + ((atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[1])].x) * (atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[1])].x));  
							newdist = newdist + ((atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[1])].y) * (atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[1])].y));  
							newdist = newdist + ((atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[1])].z) * (atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[1])].z));  
							if (newdist < dist) dist = newdist;
						
						}
						if (block[i].doublepeak == true && block[j].doublepeak == true){
							float newdist = 0;						
							newdist = newdist + ((atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[1])].x) * (atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[1])].x));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[1])].y) * (atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[1])].y));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[1])].z) * (atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[1])].z));  
							if (newdist < dist) dist = newdist;
						
						}
						if (dist > restraint){
							cout << "Found violated Restraint: " << fixed << setprecision(2) << cross[k].h1 << " " << cross[k].n1 << " -> " << cross[k].h2 << " " << cross[k].n2 << " " << setprecision(2) << cross[k].volume << endl;
							for (unsigned int l = 0; l < noes.size(); l++) if (abs(cross[k].h1 - noes[l].h1) < 0.01 && abs(cross[k].h2 - noes[l].h2) < 0.01 && abs(cross[k].n1 - noes[l].n1) < 0.1 && abs(cross[k].n2 - noes[l].n2) < 0.1) noes.erase(noes.begin() + l);
							if (optimize(atoms,noes,logfile,safety, fix) == 0){
								printnoes(atoms, noes, logfile);
								return 0;
							}
							cross.erase(cross.begin() + k);
							if (k > 0) k--;
						}
					}
					for (unsigned int k = 0; k < cross.size(); k++) if (abs(cross[k].h2 - block[i].h1) < 0.01 && abs(cross[k].n2 - block[i].n1) < 0.1) if (abs(cross[k].h1 - block[j].h2) < 0.01 && abs(cross[k].n1 - block[j].n2) < 0.1){
						float dist = 0;
						float restraint = 0;
						restraint = cross[k].dist + safety;
						restraint = restraint * restraint;
						
						dist = dist + ((atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[0])].x) * (atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[0])].x));  
						dist = dist + ((atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[0])].y) * (atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[0])].y));  
						dist = dist + ((atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[0])].z) * (atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[0])].z));  
						
						if (block[i].doublepeak == true){
							float newdist = 0;						
							newdist = newdist + ((atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[0])].x) * (atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[0])].x));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[0])].y) * (atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[0])].y));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[0])].z) * (atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[0])].z));  
							if (newdist < dist) dist = newdist;
						}
						if (block[j].doublepeak == true){
							float newdist = 0;						
							newdist = newdist + ((atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[1])].x) * (atoms[(block[i].atoms[0])].x  - atoms[(block[j].atoms[1])].x));  
							newdist = newdist + ((atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[1])].y) * (atoms[(block[i].atoms[0])].y  - atoms[(block[j].atoms[1])].y));  
							newdist = newdist + ((atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[1])].z) * (atoms[(block[i].atoms[0])].z  - atoms[(block[j].atoms[1])].z));  
							if (newdist < dist) dist = newdist;
						
						}
						if (block[i].doublepeak == true && block[j].doublepeak == true){
							float newdist = 0;						
							newdist = newdist + ((atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[1])].x) * (atoms[(block[i].atoms[1])].x  - atoms[(block[j].atoms[1])].x));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[1])].y) * (atoms[(block[i].atoms[1])].y  - atoms[(block[j].atoms[1])].y));  
							newdist = newdist + ((atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[1])].z) * (atoms[(block[i].atoms[1])].z  - atoms[(block[j].atoms[1])].z));  
							if (newdist < dist) dist = newdist;
						
						}
						if (dist > restraint){
							cout << "Found violated Restraint: " << fixed << setprecision(2) << cross[k].h1 << " " << cross[k].n1 << " -> " << cross[k].h2 << " " << cross[k].n2 << " " << setprecision(2) << cross[k].volume << endl;
							for (unsigned int l = 0; l < noes.size(); l++) if (abs(cross[k].h1 - noes[l].h1) < 0.01 && abs(cross[k].h2 - noes[l].h2) < 0.01 && abs(cross[k].n1 - noes[l].n1) < 0.1 && abs(cross[k].n2 - noes[l].n2) < 0.1) noes.erase(noes.begin() + l);
							if (optimize(atoms,noes,logfile,safety, fix) == 0){
								printnoes(atoms, noes, logfile);
								return 0;
							}
							cross.erase(cross.begin() + k);
							if (k > 0) k--;
						}
					}

				}
			}
			// Now the block is check for artefacts, if there are remaining entries in cross repeat loop
			if (cross.size() == 0) break;
		}	
		// If the last optimization was succesful best result already in logfile and artefacts are labeled
		// If not repeat loop because an artefact might be put in the block vector
		if (result == 0) return 0;
		
		cout << "Got an error in the block vector, probably caught an artefact in the block vector => ";
		cout << "Restart artefact loop and try again" << endl;
	}
	return 0;
}


unsigned int identifydp(vector<atom> & atoms, vector<noe> & noes, string logfile, float safety, vector<atom> fix){
	unsigned int maxnr  = 0;
	unsigned int maxdps = 0;
	unsigned int errorlevel = 0;

	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1){
		noes[i].crosspeaks.push_back(j);
	}



	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].crosspeaks.size() > maxnr) maxnr = noes[i].crosspeaks.size();
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) maxdps++;
	maxdps = atoms.size() - maxdps;
	/*
	while(1){
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].crosspeaks.size() == maxnr){
			cout << "Try peak: " << fixed << setprecision(2) << noes[i].h1 << " " << noes[i].n1 << " (" << fixed << setprecision(0) << noes[i].crosspeaks.size() << ")" << endl;
			noes[i].doublepeak = true;
			for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1) if (i != j) noes[j].doublepeak = true;
			errorlevel = optimize(atoms, noes, logfile, safety);	
			if (errorlevel == 0){
				 printnoes(atoms, noes, logfile);
				cout << "Identified peak as double peak. Got an errorfree optimization." << endl;
				return 0;
			}
			noes[i].doublepeak = false;
			for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1) if (i != j) noes[j].doublepeak = false;
		}
		if (maxnr == 0) break;
		maxnr--;
	}
	*/
	// Now try two doublepeaks
	for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].crosspeaks.size() > maxnr) maxnr = noes[i].crosspeaks.size();
	while(1){
		for (unsigned int i = 0; i < noes.size(); i++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (noes[i].crosspeaks.size() == maxnr){
			cout << "Try peak: " << fixed << setprecision(2) << noes[i].h1 << " " << noes[i].n1 << " (" << fixed << setprecision(0) << noes[i].crosspeaks.size() << ")" << endl;
			noes[i].doublepeak = true;
			for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1) if (i != j) noes[j].doublepeak = true;
			
			for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[i].h1 - noes[i].h2) < 0.01 && abs(noes[i].n1 - noes[i].n2) < 0.1) if (j != i){
				cout << "Try second peak: " << fixed << setprecision(2) << noes[j].h1 << " " << noes[j].n1 << " (" << fixed << setprecision(0) << noes[j].crosspeaks.size() << ")" << endl;
				noes[j].doublepeak = true;
				for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[j].h1 - noes[k].h1) < 0.01 && abs(noes[j].n1 - noes[k].n1) < 0.1) if (j != k && k != i) noes[k].doublepeak = true;
				errorlevel = optimize(atoms, noes, logfile, safety, fix);	
				if (errorlevel == 0){
					printnoes(atoms, noes, logfile);
					cout << "Identified peaks as double peaks. Got an errorfree optimization." << endl;
					return 0;
				}
				noes[j].doublepeak = false;
				for (unsigned int k = 0; k < noes.size(); k++) if (abs(noes[j].h1 - noes[k].h1) < 0.01 && abs(noes[j].n1 - noes[k].n1) < 0.1) if (j != k && k != i) noes[k].doublepeak = false;
				
			}
			noes[i].doublepeak = false;
			for (unsigned int j = 0; j < noes.size(); j++) if (abs(noes[i].h1 - noes[j].h1) < 0.01 && abs(noes[i].n1 - noes[j].n1) < 0.1) if (i != j) noes[j].doublepeak = false;
		}
		if (maxnr == 0) break;
		maxnr--;
	}
	
	
	
	
	
	
	
	
	return 0;
}							  
							  
							  
							  
							  
							  
							  
							  
							  
							  
							  
