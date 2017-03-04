#pragma once
#ifndef _FILEOPERATION_H
#define _FILEOPERATION_H

#include "../include/LDPC_Coding.h"
#define INFILE_NAME "data/datain.txt"
#define OUTFILE_NAME "data/dataout.txt"

/*API for reading data from the file*/
void OpenDataFile();
void CloseDataFile();
int ReadData(int *buf);  //return total number of data in the file
int GetIdx(int *buf); //return the index data from the buffer at the current cursor
int ReadLineData(int *buf); //return the total number of data in a sigle line, 
							//the buffer stores the index data from a single line

void WriteMatrixData(int *buf, int row, int col);   //write the matrix data to the txt file
//int WriteVectorData(LDPC_int *buf, int length);

void WriteDecodedData(LDPC_int *buf, int length);


#endif