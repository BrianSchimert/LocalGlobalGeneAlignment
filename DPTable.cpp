
#include <iostream>
#include "DPTable.h"
#include <vector>
#include <string>
#include <stack>
using namespace std;



vector<vector<DPCell>> T;
int m_max, n_max, match_score, mismatch_score, h_score, g_score;
string seq1, seq2;
int numMatches=0, numMismatch=0, numGaps=0;
stack<char> s1ans;
stack<char> s2ans;



void DPTable::initGlobal(int m, int n, int match, int mismatch, int h, int g, string s1, string s2)
{
	seq1 = s1;
	seq2 = s2;
	m_max = m;
	n_max = n;
	match_score = match;
	mismatch_score = mismatch;
	h_score = h;
	g_score = g;
	for (int i = 0; i <= m ; i++) {
		vector<DPCell> row;
		for (int j = 0; j <= n; j++) {
			DPCell cell;
			//cell.dScore = 0;
			//cell.iScore = 0;
			//cell.sScore = 0;
			cell.score = 0;
			row.push_back(cell);
		}
		T.push_back(row);
	}
	T[0][0].sScore = 0;
	T[0][0].iScore = 0;
	T[0][0].dScore = 0;
	for (int i = 1; i <= m ; i++) {
		T[i][0].iScore = -1000;
		T[i][0].dScore = h + i*g;
		T[i][0].sScore = -1000;
		T[i][0].score = i * g;
		//T[i][0].direction_case = "deletion";// or delition?
	}
	for (int j = 1; j <= n; j++) {
		T[0][j].iScore = h + j*g;
		T[0][j].dScore = -1000;
		T[0][j].sScore = -1000;
		T[0][j].score = j * g;
		//T[0][j].direction_case = "insertion";   //or insertion?
	}
	
}

void DPTable::initLocal(int m, int n, int match, int mismatch, int h, int g, string s1, string s2)
{
	seq1 = s1;
	seq2 = s2;
	m_max = m;
	n_max = n;
	match_score = match;
	mismatch_score = mismatch;
	h_score = h;
	g_score = g;
	for (int i = 0; i <= m; i++) {
		vector<DPCell> row;
		for (int j = 0; j <= n; j++) {
			DPCell cell;
			//cell.dScore = 0;
			//cell.iScore = 0;
			//cell.sScore = 0;
			cell.score = 0;
			row.push_back(cell);
		}
		T.push_back(row);
	}
	T[0][0].sScore = 0;
	T[0][0].iScore = 0;
	T[0][0].dScore = 0;
	for (int i = 1; i <= m; i++) {
		T[i][0].iScore = -1000;
		T[i][0].dScore = 0;
		T[i][0].sScore = -1000;
		T[i][0].score = 0;
		//T[i][0].direction_case = "deletion";
	}
	for (int j = 1; j <= n; j++) {
		T[0][j].iScore = 0;
		T[0][j].dScore = -1000;
		T[0][j].sScore = -1000;
		T[0][j].score = 0;
		//T[0][j].direction_case = "insertion";  
	}

}



void DPTable::fillGlobalTable() {
	int try1;
	int try2;
	int try3;
	for (int i = 1; i <= m_max; i++) {
		for (int j = 1; j <= n_max; j++) {
			//Starting to find and set s score for current cell
			if (seq1[i] == seq2[j]) {  //match
				try1 = T[i - 1][j - 1].iScore + match_score;
				try2 = T[i - 1][j - 1].dScore + match_score;
				try3 = T[i - 1][j - 1].sScore + match_score;
			}
			else {	//mismatch
				try1 = T[i - 1][j - 1].iScore + mismatch_score;
				try2 = T[i - 1][j - 1].dScore + mismatch_score;
				try3 = T[i - 1][j - 1].sScore + mismatch_score;
			}
			// set S score for current cell based on diagonal predecessors' max value (I, D, S)
			T[i][j].sScore = maxOf(try1, try2, try3);  

			//Starting to find and set D score for current cell
			try1 = T[i - 1][j].dScore + g_score;
			try2 = T[i - 1][j].sScore + g_score + h_score;  //Deletion
			try3 = T[i - 1][j].iScore + g_score + h_score;

			T[i][j].dScore = maxOf(try1, try2, try3);

			//Starting to find and set I scores for current cell
			try1 = T[i][j - 1].iScore + g_score;
			try2 = T[i][j - 1].sScore + h_score + g_score;
			try3 = T[i][j - 1].dScore + h_score + g_score;

			T[i][j].iScore = maxOf(try1, try2, try3);
		}
	}
}


void DPTable::fillLocalTable() {
	int try1;
	int try2;
	int try3;
	int checkGTZero;
	for (int i = 1; i <= m_max; i++) {
		for (int j = 1; j <= n_max; j++) {
			//Starting to find and set s score for current cell
			if (seq1[i] == seq2[j]) {  //match
				try1 = T[i - 1][j - 1].iScore + match_score;
				try2 = T[i - 1][j - 1].dScore + match_score;
				try3 = T[i - 1][j - 1].sScore + match_score;
			}
			else {	//mismatch
				try1 = T[i - 1][j - 1].iScore + mismatch_score;
				try2 = T[i - 1][j - 1].dScore + mismatch_score;
				try3 = T[i - 1][j - 1].sScore + mismatch_score;
			}
			// set S score for current cell based on diagonal predecessors' max value (I, D, S)
			checkGTZero = maxOf(try1, try2, try3);
			if (checkGTZero > 0) {
				T[i][j].sScore = checkGTZero;
			}
			else {
				T[i][j].sScore = 0;
			}
			//Starting to find and set D score for current cell
			try1 = T[i - 1][j].dScore + g_score;
			try2 = T[i - 1][j].sScore + g_score + h_score;  //Deletion
			try3 = T[i - 1][j].iScore + g_score + h_score;

			checkGTZero = maxOf(try1, try2, try3);
			if (checkGTZero > 0) {
				T[i][j].dScore = checkGTZero;
			}
			else {
				T[i][j].dScore = 0;
			}

			//Starting to find and set I scores for current cell
			try1 = T[i][j - 1].iScore + g_score;
			try2 = T[i][j - 1].sScore + h_score + g_score;
			try3 = T[i][j - 1].dScore + h_score + g_score;

			checkGTZero = maxOf(try1, try2, try3);
			if (checkGTZero > 0) {
				T[i][j].iScore = checkGTZero;
			}
			else {
				T[i][j].iScore = 0;
			}
		}
	}
}

void DPTable::retrace(int x, int y) {
	//int i = m_max; 
	//int j = n_max;
	int i = x;
	int j = y;
	DPCell currCell = T[i][j];
	int maxCellScore;
	while (i >= 0 && j >= 0) {
		currCell = T[i][j];
		maxCellScore = maxOf(currCell.iScore, currCell.dScore, currCell.sScore);
		if (maxCellScore == currCell.sScore) { //substitution case, go diagonal (and check match/mismatch)
		   // add ai bj to path, i-- j--;
			if (seq1[i] == seq2[j]) {
				numMatches++;
			}
			else {
				numMismatch++;
			}
			s1ans.push(seq1[i]);
			s2ans.push(seq2[j]);
			i--;
			j--;
		}
		else if (maxCellScore == currCell.iScore) { //insertion case, go left
			//add insertion gap to path - bj, j--;
			s1ans.push('_');
			s2ans.push(seq2[j]);
			j--;
			numGaps++;
		}
		else if (maxCellScore == currCell.dScore) { //deletion case, go up
		   //add deletion gap to path ai - , i--;
			s1ans.push(seq1[i]);
			s2ans.push('_');
			i--;
			numGaps++;
		}

	}

	cout << seq1 << "\n";
	cout << seq2 << "\n";
	cout << "\n\n";
	int s1size = s1ans.size();
	for (int x = 0; x < s1size; x++) {
		cout << s1ans.top();
		s1ans.pop();
	}
	cout << "\n\n";
	int s2size = s2ans.size();
	for (int x = 0; x < s2size; x++) {
		cout << s2ans.top();
		s2ans.pop();
	}

	cout << "\n\n\n\n";
	cout << "Number of matches: " << numMatches << " Mismatches: " << numMismatch << "  Gaps: " << numGaps << "\n\n";
	//cout << seq1[0];
	//cout << seq2[0];
}


int DPTable::maxOf(int x, int y, int z) {
	if (x >= y && x >= z) {
		return x;
	}
	else if (y >= x && y >= z) {
		return y;
	}
	else if (z >= x && z >= y) {
		return z;
	}
}


void DPTable::printBestScore(int i, int j) {
	DPCell answer = T[i][j];
	int maxScore = maxOf(answer.dScore, answer.iScore, answer.sScore);
	cout << maxScore;

}

void DPTable::printTable() {
	int counter = 0;
	int max;
	for (int i = 0; i < m_max; i++) {
		
		if (counter % n_max == 0) {
			cout << "\n\n";
		}
		cout << "{";
		for (int j = 0; j < n_max; j++) {
			max = maxOf(T[i][j].sScore, T[i][j].dScore, T[i][j].iScore);
			cout << max << " ";
			counter++;
		}
		cout << "}";
	}

}


//
//void DPTable::printTable() {
//	int max;
//	int counter = 0;
//	for (int i = 0; i < T.size(); i++)
//	{
//		if (counter % n_max == 0) {
//			cout << "\n\n";
//		}
//		for (int j = 0; j < T[i].size(); j++)
//		{
//			max = maxOf(T[i][j].sScore, T[i][j].dScore, T[i][j].iScore);
//			cout << max << " ";
//		}
//	}
//}




//alternative init with dynamically allocated 2d array, instead of vector of vectors

//void DPTable::init(int m, int n, int match, int mismatch, int h, int g)
//{
	//cout << m;
	//cout << n;
	//DPCell[m][n] dpTable;
	//int** T = new int* [m];
	//for (int i = 0; i < m; ++i) {
	//	T[i] = new int[n];
	//}
	//
	//T[0][0] = 0;
	//for (int i = 0; i < m; i++) {
	//	T[i][0] = i * g;
	//}
	//for (int j = 0; j < n; j++) {
	//	T[0][j] = j * g;
	//}
	//for (int i = 0; i < m; i++) {
	//	for (int j = 0; j < n; j++) {
	//		cout << T[i][j];
	//	}
	//}
//}


//alterntive fillTable function with regular linear gap model
//void DPTable::fillTable() {
//	int try1;
//	int try2;
//	int try3;
//	for (int i = 1; i < m_max; i++) {
//		for (int j = 1; j < n_max; j++) {
//			if (seq1[i] == seq2[j]) {  //match
//				try1 = T[i - 1][j - 1].iScore + match_score;
//				try2 = T[i - 1][j - 1].dScore + match_score;
//				try3 = T[i - 1][j - 1].sScore + match_score;
//			}
//			else {	//mismatch
//				try1 = T[i - 1][j - 1].iScore - mismatch_score;
//				try2 = T[i - 1][j - 1].dScore - match_score;
//				try3 = T[i - 1][j - 1].sScore - match_score;
//			}
//
//			try2 = T[i - 1][j].score + g_score;
//			try3 = T[i][j - 1].score + g_score;
//			if (try1 >= try2 && try1 >= try3) {
//				T[i][j].score = try1;
//				T[i][j].direction_case = "substitution";
//			}
//			else if (try2 >= try1 && try2 >= try3) {
//				T[i][j].score = try2;
//				T[i][j].direction_case = "deletion"; //  ???????????
//			}
//			else if (try3 >= try2 && try3 >= try1) {
//				T[i][j].score = try3;
//				T[i][j].direction_case = "insertion";  //  ????? 
//			}
//		}
//	}
//}