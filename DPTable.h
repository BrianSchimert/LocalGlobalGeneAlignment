#pragma once
#include "DPCell.h"
#include <string>

	class DPTable {

	public:
		void initGlobal(int m, int n, int match, int mismatch, int h, int g, std::string s1, std::string s2);
		void initLocal(int m, int n, int match, int mismatch, int h, int g, std::string s1, std::string s2);
		//void update();
		void printTable();
		void fillGlobalTable();
		void fillLocalTable();
		int maxOf(int x, int y, int z);
		void printBestScore(int i, int j);
		void retrace(int x, int y);
	};


