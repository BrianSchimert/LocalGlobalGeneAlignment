#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "GeneSequence.h"
#include "DPTable.h"


using namespace std;

vector<GeneSequence> genes;

void openParameterFile(string myfile, int &match, int &mismatch, int &h, int &g) {
	ifstream myFile(myfile);
	int index = 0;
	string line;
	while (getline(myFile, line)) {
		istringstream streamer(line);
		string word;
		vector<string> thisline;
		while (streamer >> word) {
			thisline.push_back(word);
		}
		if (index == 0) {
			match = stoi(thisline[1]);
		}
		else if (index == 1) {
			mismatch = stoi(thisline[1]);
		}
		else if (index == 2) {
			h = stoi(thisline[1]);
		}
		else if (index == 3) {
			g = stoi(thisline[1]);
		}
		index++;
	}
}

void openGenomeFile(string myfile) {
	ifstream myFile(myfile);
	string line;
	//GeneSequence currentGene;
	GeneSequence currentGene;
	int time = 0;
	while (getline(myFile, line)) {
		if (line.empty()) {
			continue;
		}
		else if (line[0] == '>') { //new sequence found
			if (time != 0) {
				genes.push_back(currentGene);
				currentGene.sequence.clear();
			}	
			string name = line.substr(1, line.find_first_of(" "));
			currentGene.name = name;
			time = 1;
		}
		else {  // (sequence lines) ie. lines that are not part of the name, or empty
			currentGene.sequence += line;
		}

	}
	genes.push_back(currentGene);

	//for (int i = 0; i < genes.size(); i++) {
	//	cout << genes.at(i).sequence << "\n";
	//}

}

// <executable name> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>
int main(int argc, char* argv[])
{

	ofstream cout("GlobalAlignment.txt");
	int match, mismatch, h, g;
	int local_flag = 0;
	string genomeFile;
	if (argc > 3) {
		//for (int i = 0; i < argc; i++)
		//	std::cout << argv[i] << "\n";
		genomeFile = argv[1];
		local_flag = stoi(argv[2]);
		string paramFile = argv[3];
		openParameterFile(paramFile, match, mismatch, h, g);
		openGenomeFile(genomeFile);
		
	}
	else {
		genomeFile = "input.fasta";
		//genomeFile = "Opsin1_colorblindness_gene.fasta";
		//local_flag = 0;  //global alignment by default
		openParameterFile("parameters.config", match, mismatch, h, g);
		openGenomeFile(genomeFile);

		}
	string s1 = genes.at(0).sequence;
	string s2 = genes.at(1).sequence;
	string s1Name = genes.at(0).name;
	string s2Name = genes.at(1).name;
	int m = s1.length();
	int n = s2.length();
	
	
	cout << "Sequence 1 Info: " << "\n";
	cout << "Name: " << s1Name << "   Length: " << m;
	cout << "\n\n";
	//cout << s1;
	cout << "\n\n";

	cout << "Sequence 2 Info: " << "\n";
	cout << "Name: " << s2Name << "   Length: " << n;
	cout << "\n\n";
	//cout << s2;
	cout << "\n\n";
	DPTable dpTable;
	
	if (local_flag == 0) {
		dpTable.initGlobal(m, n, match, mismatch, h, g, s1, s2);
		dpTable.fillGlobalTable();
		
	}
	else {
		dpTable.initLocal(m, n, match, mismatch, h, g, s1, s2);
		dpTable.fillLocalTable();
	}

	cout << "Parameters: " << "Match: " << match << "   Mismatch: " << mismatch << "   h: " << h << "   g: " << g << "\n\n";
	if (local_flag == 0) {
		cout << "Global Optimal Score:  ";
		dpTable.printBestScore(m, n);
		cout << "\n\n";
		dpTable.retrace(m, n);
	}
	else {
		cout << "Local Optimal Score: ";
		dpTable.printBestScore(m, n);
		cout << "\n\n";
		dpTable.retrace(m,n);
	}

	dpTable.printTable();
	

	


	return 0;
}
