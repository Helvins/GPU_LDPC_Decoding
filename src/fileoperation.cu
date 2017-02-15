#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include "../include/fileoperation.h"
using namespace std;


static int cursor = -1;
ifstream infile;
ofstream outfile;
string line;

void OpenDataFile() {
	infile.open(INFILE_NAME, ios::in);
	if (!infile.is_open()) {
		cout << "Error opening infile!" << endl;
		exit(1);
	}
	outfile.open(OUTFILE_NAME, ios::out | ios::trunc);      //open output file will cover the former one
	if (!outfile.is_open()) {
		cout << "Error opening outfile!" << endl;
		exit(1);
	}
}

void CloseDataFile() {
	infile.close();
	outfile.close();
}

int ReadData(int *buf) {
	memset(buf, -1, 15); //initialize the memory space
	int digitcount = 0, sum = 0;
	if (infile.is_open()) {
		int totalcount = 0;
		while (!infile.eof()) {
			getline(infile, line);
			do {
				if (line[digitcount] >= '0' && line[digitcount] <= '9') {
					//cout << line[digitcount] << " ";
					sum = sum * 10 + (line[digitcount] - '0');
				}

				else {
					//condition of space	
					buf[totalcount] = sum;
					totalcount++;
					sum = 0;
				}

				digitcount++;
				
				if (line[digitcount] == '\0') {
					//the last number in the line
					buf[totalcount] = sum;
					totalcount++;
				}

			} while (line[digitcount] != '\0');
			/*
			string::const_iterator p, q = line.begin(), end = line.end();
			while ((p = find_if(q, end, isdigit)) != end) {
				q = find_if_not(p, end, isdigit);
				buf[totalcount] = stoul(line.substr(distance(line.cbegin(), p), distance(p, q)));
				//cout << buf[totalcount] << " ";
				totalcount++;
			}*/
		}
		return totalcount;
	}
	else {
		cout << "Error opening file!" << endl;
		return -1;
	}
}

int GetIdx(int *buf) {
	cursor++;
	//cout << "Current cusor is: " << cursor << " ";
	return buf[cursor];
}

int ReadLineData(int *buf) {
	memset(buf, -1, 15); //initialize the memory space
	int totalcount = 0, digitcount = 0, sum = 0;
	if (!infile.eof()) {
		getline(infile, line);
		do {
			if(line[digitcount]>='0' && line[digitcount]<='9') {
				sum  = sum*10+(line[digitcount]-'0');	
			}
			
			else {
				buf[totalcount] = sum;
				totalcount++;
				sum = 0;
			}
			
			digitcount++;
			/*
			if (line[digitcount] == '\0') {
				//the last number in the line
				buf[totalcount] = sum;
				totalcount++;
			}*/

		} while (line[digitcount] != '\0');

		/*
		getline(infile, line);

		string::const_iterator p, q = line.begin(), end = line.end();
		while ((p = find_if(q, end, isdigit)) != end) {
			q = find_if_not(p, end, isdigit);
			buf[totalcount] = stoul(line.substr(distance(line.cbegin(), p), distance(p, q)));
			totalcount++;
		}*/
		return totalcount;
	}
	else {
		cout << "No extra lines for reading in the file!" << endl;
		return -1;
	}
	
}

void WriteMatrixData(int *buf, int row, int col) {
	outfile << "Now writing new matrix data"<<endl;
	for (int i = 0;i < row;i++) {
		for (int j = 0;j < col; j++) {
			outfile << buf[i*col+j] << " ";
		}
		outfile << endl;
	}
	outfile<<endl<<endl;
}



