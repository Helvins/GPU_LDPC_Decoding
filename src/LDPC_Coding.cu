/* /O2 Compiling optimized removing /Zl->Zi /RTC option */

#include "../include/LDPC_Coding.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
#include <cmath>
#include <stddef.h>
#include "../include/fileoperation.h"
using namespace std;


int idxdata[MAX_COLUMN_WEIGHT];
/*the order is CodeWord_Size, Info_Size, d_pitch_r, d_pitch_c*/
__constant__ int CUDA_PARA_SET[4];
cudaError_t ret;

/*Constructor without passing parameters*/
LDPC_Coding::LDPC_Coding() {
	CodeWord_Size = DEFAULT_CODEWORD_SIZE;
	Code_Rate = DEFAULT_CODE_RATE;
	Info_Size = (unsigned int)(DEFAULT_CODEWORD_SIZE*DEFAULT_CODE_RATE);
	Blk_Num = DEFAULT_BLOCK_NUM;
}

/*Overloading the constructor*/
LDPC_Coding::LDPC_Coding(int size, float rate) {
	CodeWord_Size = size;
	Code_Rate = rate;
	Info_Size = (unsigned int)(size*rate);
	Blk_Num = (unsigned int)(Info_Size / 360.0);

}

bool LDPC_Coding::Memory_Space_Allocation() {
	/*allocate memory space for Info_Size*CodeWord_Size check matrix*/
	Check_Matrix = (LDPC_int *)malloc(Info_Size*CodeWord_Size*sizeof(LDPC_int));				//warning:can be substituted by cudaHostMalloc
	if (!Check_Matrix) {
		printf("Error:Can't allocate memory space for check matrix!\n");
		return FAIL;
	}

	
	/*use page-locked memory in the host*/
	ret = cudaMallocHost((void**)&Index_Row_Matrix, Info_Size*MAX_ROW_WEIGHT*sizeof(int));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		printf("Error:Can't allocate memory space for row index matrix!\n");
		return FAIL;
	}

	/* allocate memory space for column index matrix (MAX_COLUMN_WEIGHT*CodeWord_Size) */
	/*
	Index_Col_Matrix = (int *)malloc(MAX_COLUMN_WEIGHT*CodeWord_Size*sizeof(int));
	if (!Index_Col_Matrix) {
			printf("Error:Can't allocate memory space for column index matrix!\n");
			return FAIL;
	}*/
	/*use page-locked memory in the host*/
	ret = cudaMallocHost((void**)&Index_Col_Matrix, MAX_COLUMN_WEIGHT*CodeWord_Size*sizeof(int));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		printf("Error:Can't allocate memory space for row index matrix!\n");
		return FAIL;
	}

	/*initialization for a K*N check matrix */
	for (int row = 0;row < Info_Size;row++) {
		for (int col = 0;col < CodeWord_Size; col++) {
			//Check_Matrix[row][col] = LDPC_0;
			Check_Matrix[row*CodeWord_Size+col] = LDPC_0;
		}
	}

	/*initialization for the row index matrix*/
	for (int h = 0;h < Info_Size;h++) {
		for (int i = 0;i < MAX_ROW_WEIGHT;i++) {
			//Index_Row_Matrix[h][i] = -1;
			Index_Row_Matrix[h*MAX_ROW_WEIGHT+i] = -1;
		}
	}

	/*initialization for the column index matrix*/
	for (int h = 0;h < MAX_COLUMN_WEIGHT;h++) {
		for (int i = 0;i < CodeWord_Size;i++) {
			//Index_Col_Matrix[h][i] = -1;
			Index_Col_Matrix[h*CodeWord_Size+i] = -1;
		}
	}

	return SUCCESS;
}

/*Destructor*/
LDPC_Coding::~LDPC_Coding() {
	/*recollect the memory space in the Check Matrix*/
	//free(Check_Matrix);

	/*recollect the memory space in the credibility array p0 and p1*/
	/*
	free(p0);
	free(p1);
	*/
	/*recollect the momory space in the Row Index Matrix*/
	/*
	free(Index_Row_Matrix);
	free(r0);
	free(r1);
	*/
	/*recollect the momory space in the Column Index Matrix*/
	/*
	free(Index_Col_Matrix);
	free(q0);
	free(q1);
	*/
	/*close data files including input and output txt files*/
	CloseDataFile(); 
}

bool LDPC_Coding::Sparse_Cyclic_Matrix_Construct() {
	register unsigned short trrow, trcol;
	int totalcount = 0, currentcount = 0;
	OpenDataFile();

	for (short subblock_num = 0;subblock_num < Blk_Num; subblock_num++) {
		totalcount = ReadLineData(idxdata);  // the index data is stored in the idxdata buffer	
		if (totalcount != -1) {
			for (short col = 0; col < 360; col++) {
				currentcount = 0;

				while (currentcount < totalcount) {
					//-1 stands for the EOF
	
					trrow = (idxdata[currentcount] + Blk_Num * col) % (CodeWord_Size - Info_Size);
					trcol = col + subblock_num * 360;

					//Check_Matrix[trrow][trcol] = LDPC_1;
					Check_Matrix[trrow*CodeWord_Size+trcol] = LDPC_1;
					//Index_Col_Matrix[currentcount][trcol] = trrow;   //the 2-dimensional coordinate of non-zero element is (trrow,trcol) 
					Index_Col_Matrix[currentcount*CodeWord_Size+trcol] = trrow; 
					currentcount++;
				}
				//system("pause");

			}
		}
		else if (totalcount == -1 && subblock_num < Blk_Num) {
			cout << "Errot reading file data" << endl;
			return FAIL;
		}
	}

	return SUCCESS;
}

/*Construct the double diagonal matrix in the check matrix*/
bool LDPC_Coding::Double_Diagonal_Matrix_Construct() {
	register int mainrow, maincol, auxrow, auxcol;
		for (int diag = 0;diag < (CodeWord_Size - Info_Size);diag++) {
			mainrow = diag;
			maincol = Info_Size + diag;

			auxrow = mainrow + 1;
			auxcol = maincol;

			if (diag != (CodeWord_Size - Info_Size - 1)) {	
				Check_Matrix[mainrow*CodeWord_Size+maincol] = LDPC_1;
				Check_Matrix[auxrow*CodeWord_Size+auxcol] = LDPC_1;
				//Check_Matrix[mainrow][maincol] = LDPC_1;
				//Check_Matrix[auxrow][auxcol] = LDPC_1;
				
				/*update the column index matrix simultaneously*/
				Index_Col_Matrix[maincol] = mainrow;
				Index_Col_Matrix[CodeWord_Size+auxcol] = auxrow;
				//Index_Col_Matrix[0][maincol] = mainrow;
				//Index_Col_Matrix[1][auxcol] = auxrow;
			}
			else {										  //the auxiliary row should change to the first row when the diag equals to the last column 
				Check_Matrix[mainrow*CodeWord_Size+maincol] = LDPC_1;  //the auxiliary diagonal is out of bound
				//Check_Matrix[mainrow][maincol] = LDPC_1;  //the auxiliary diagonal is out of bound

				//Check_Matrix[0][CodeWord_Size - 1] = LDPC_1;   //the element on the top right corner is 1
				Index_Col_Matrix[maincol] = mainrow;
				//Index_Col_Matrix[0][maincol] = mainrow;
				//Index_Col_Matrix[1][auxcol] = 0;
			}
		}
		
		return SUCCESS;
}

bool LDPC_Coding::Check_Matrix_Construct() {						// the size of check matrix is K*N
	if (Sparse_Cyclic_Matrix_Construct() && Double_Diagonal_Matrix_Construct()) {
		printf("Constuct check matrix success\n");
		return SUCCESS;
	}
	else
		return FAIL;
}   

bool LDPC_Coding::Calculate_Max_Row_Weight() {
	int max = 0, tmp = 0;
	for (int i = 0;i < Info_Size;i++) {
		tmp = 0;
		for (int j = 0;j < CodeWord_Size;j++) {
			//if (Check_Matrix[i][j]) {
			if (Check_Matrix[i*CodeWord_Size+j]) {
				tmp++;
			}
		}
		max = (tmp > max) ? tmp : max;
	}
	cout << "The max row weight is " << max <<endl;
	return SUCCESS;
}

/* the memory space of check matrix will be freed after this funcion executed */
bool LDPC_Coding::Scan_Check_Matrix() {
	
	int currentcount;
	for (int i = 0;i < Info_Size;i++) {
		currentcount = 0;
		for (int j = 0;j < CodeWord_Size;j++) {
			//if (Check_Matrix[i][j]) {
			if (Check_Matrix[i*CodeWord_Size+j]) {
				Index_Row_Matrix[i*MAX_ROW_WEIGHT+currentcount] = j;                   //record the column position
				//Index_Row_Matrix[i][currentcount] = j;                  
				currentcount++;
			}
		}
	}

	/*recollect the memory space in the Check Matrix*/
	free(Check_Matrix);
	/*
	for (int i = 0;i < Info_Size;i++) {
		delete[]Check_Matrix[i];
		//free((void*)Check_Matrix[i]);
	}
	delete[]Check_Matrix;
	*/
	
	return SUCCESS;
}

void LDPC_Coding::WriteData(){
	printf("Now try to write the row and column index matrix into the output file\n");
	WriteMatrixData(Index_Row_Matrix, Info_Size, MAX_ROW_WEIGHT);
	WriteMatrixData(Index_Col_Matrix, MAX_ROW_WEIGHT, CodeWord_Size);
}

/*LDPC encoding to generate a code word*/
bool LDPC_Coding::LDPC_Encoding(LDPC_int *info_seq, LDPC_int *code_word) {
	register int currentcount, tmp;
	/*step 1 construct the parity check part*/
	for (int num = 0;num < Info_Size;num++) {
		code_word[num] = info_seq[num];										//the code word has the same information sequence from the format K bits
		currentcount = 0;
		//while (Index_Col_Matrix[currentcount][num] != -1) {
		while (Index_Col_Matrix[currentcount*CodeWord_Size+num] != -1) {
			tmp = Index_Col_Matrix[currentcount*CodeWord_Size+num]+Info_Size; 
			//tmp = Index_Col_Matrix[currentcount][num]+Info_Size;  		//calculate the absolute position of the non-zero element in the num th column		
			code_word[tmp] ^= info_seq[num];								//process the xor operation to update the parity bit based on the information bit      
			currentcount++;
		}
	}
	for (int num = Info_Size+1;num < CodeWord_Size;num++) {					//process the xor operation between neighbouring parity bits
		code_word[num] ^= code_word[num - 1];
	}
	//code_word[Info_Size] ^= code_word[CodeWord_Size - 1];
	return SUCCESS;
}

bool LDPC_Coding::LDPC_Encoding_Check(LDPC_int *code_word) {
	LDPC_int sum;
	int count = 0;
	LDPC_int *code_word_cp = new LDPC_int[CodeWord_Size];
	
	/* using the whole check matrix */
	/*
#pragma ompparallel for 
	for (int i = 0;i < Info_Size;i++) {                    //every row can be processed in parallel
		sum = false;
		for (int j = 0;j < CodeWord_Size;j++) {
			sum ^= (Check_Matrix[i][j] ^ code_word[j]);
		}
		if (sum) {
			cout << "LDPC Encoded Check fail" << endl;
			return FAIL;
		}
	}
		cout << "LDPC Encoded Check success" << endl;
		return SUCCESS;
		*/
	
	/*using the index matrix*/
	int number = 0;
	for (int i = 0;i < Info_Size;i++) {
		memcpy(code_word_cp, code_word, CodeWord_Size * sizeof(LDPC_int));
		sum = LDPC_0;
		count = 0;
		while (Index_Row_Matrix[MAX_ROW_WEIGHT*i+count] != -1) {                        //pick out the non-zero element
		//while (Index_Row_Matrix[i][count] != -1) { 
			code_word_cp[Index_Row_Matrix[MAX_ROW_WEIGHT*i+count]] ^= LDPC_1;
			//code_word_cp[Index_Row_Matrix[i][count]] ^= LDPC_1;
			count++;
		}


		for (int j = 0;j < CodeWord_Size;j++) {
			sum ^= code_word_cp[j];
		}

		if (sum) {
			number++;
			return FAIL;
		}
		
	}
	printf("The error bit number is: %d\n", number);
		return SUCCESS;
}


float LDPC_Coding::ErrorRate_Check(LDPC_int *info_seq, LDPC_int *de_info_seq) {
	int sum = 0;
	double error_rate;
	for(int i =0;i<Info_Size;i++){
		if (info_seq[i] != de_info_seq[i]) {
			//printf("%d ", i);
			sum++;
		}
	}
	error_rate = (float)sum / Info_Size;
	return error_rate;
}

void LDPC_Coding::Rand_Seq_Generator(LDPC_int *info_seq) {
	//int temp;
	/*for (int i = 0;i < Info_Size;i++) {
		temp = rand() % 2;  // randomly generate an integer 0 or 1
		info_seq[i] = (temp == 0) ? LDPC_0 : LDPC_1;
	}*/

	for (int i = 0;i < Info_Size/2;i++) {
		
		info_seq[i] = LDPC_1;
	}
	for (int i = Info_Size/2;i < Info_Size;i++) {
		info_seq[i] = LDPC_0;
	}

}



/*implementation of functions in class LDPC_Coding_d*/
LDPC_Coding_d::LDPC_Coding_d(){
	CodeWord_Size = DEFAULT_CODEWORD_SIZE;
	Code_Rate = DEFAULT_CODE_RATE;
	Info_Size = (unsigned int)(DEFAULT_CODEWORD_SIZE*DEFAULT_CODE_RATE);
	Col_Weight = MAX_COLUMN_WEIGHT;
	Row_Weight = MAX_ROW_WEIGHT;
}

LDPC_Coding_d::~LDPC_Coding_d(){
	
}

bool LDPC_Coding_d::Devide_Memory_Space_Allocation(){
	/*allocating memory space for waveform array*/
	ret = cudaMalloc((void**)&waveform_d, CodeWord_Size*sizeof(float));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	/*allocating dynamic memory space for array p0_d,p1_d,p0_init_d and p1_init_d*/
	ret = cudaMalloc((void**)&p0_d, CodeWord_Size*sizeof(float));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}/*
	else
		printf("The device pointer p0_d is %p\n", p0_d);*/
	ret = cudaMalloc((void**)&p1_d, CodeWord_Size*sizeof(float));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	ret = cudaMalloc((void**)&p0_init_d, CodeWord_Size*sizeof(float));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	ret = cudaMalloc((void**)&p1_init_d, CodeWord_Size*sizeof(float));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	
	/*allocating dynamic memory space for array decoded sequence*/
	ret = cudaMalloc((void**)&de_info_seq_d, CodeWord_Size*sizeof(LDPC_int));
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	
	/*allocating dynamic memory space for the row and column index mapping array*/
	/*while using cudaMallocPicth for allocating memory, the actual address can be calculated as follow(2D array -> 1D linear array*/
	/* T* pElement = (T*)((char*)BaseAddress+Row*pitch)+column */
	ret = cudaMallocPitch((void**)&Index_Row_Matrix_d, &d_pitch_r, Row_Weight*sizeof(int), Info_Size);
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	else{
		printf("The row pitch in the device is: %d\n", int(d_pitch_r));
	}
	
	ret = cudaMallocPitch((void**)&Index_Col_Matrix_d, &d_pitch_c, CodeWord_Size*sizeof(int), Col_Weight);
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	else{
		printf("The column pitch in the device is: %d\n", int(d_pitch_c));
	}
	
	/*allocating dynamic memory space for the forward and backward probability matrix r0,r1,q0 and q1*/
	ret = cudaMallocPitch((void**)&r0_d, &d_pitch_r, Row_Weight*sizeof(float), Info_Size);
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	
	ret = cudaMallocPitch((void**)&r1_d, &d_pitch_r, Row_Weight*sizeof(float), Info_Size);
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}

	ret = cudaMallocPitch((void**)&q0_d, &d_pitch_c, CodeWord_Size*sizeof(float), Col_Weight);
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}
	
	ret = cudaMallocPitch((void**)&q1_d, &d_pitch_c, CodeWord_Size*sizeof(float), Col_Weight);
	checkCudaErrors(ret);
	if(ret != cudaSuccess){
		//printf("Error allocating memory!\n");
		return FAIL;
	}

	return SUCCESS;
}

bool LDPC_Coding_d::CUDA_Configuration(dim3 cfg_para[]){
	/*create a grid containing 256 blocks with 256 threads of each block*/
	dim3 numBlocksP1(DEFAULT_CUDA_BLOCK_NUM*4,1,1);
	dim3 threadsPerBlockP1(MAX_ROW_WEIGHT, 32, 1);
	/*assign the size of shared memory in kernel function P1*/
	dim3 size_memP1(1, 1, 1);
	
	dim3 numBlocksP2(DEFAULT_CUDA_BLOCK_NUM*4,1,1);
	dim3 threadsPerBlockP2(MAX_COLUMN_WEIGHT, 64, 1);
	dim3 size_memP2(1024, 1, 1);
	
	cfg_para[0] = numBlocksP1;
	cfg_para[1] = threadsPerBlockP1;
	cfg_para[2] = size_memP1;
	
	cfg_para[3] = numBlocksP2;
	cfg_para[4] = threadsPerBlockP2;
	cfg_para[5] = size_memP2;
	
	return SUCCESS;
}


bool LDPC_Coding_d::CUDA_Memcpy_todev(LDPC_Coding entity, bool constant, float *waveform){
	//prinf("Index_Col_Matrix: %d,")
	if(constant){																//only transfer Index matrix 
		/*copy index matrix data*/
		ret = cudaMemcpy2D(Index_Row_Matrix_d, d_pitch_r, entity.Index_Row_Matrix, sizeof(int)*MAX_ROW_WEIGHT, sizeof(int)*MAX_ROW_WEIGHT, Info_Size, cudaMemcpyHostToDevice);
		checkCudaErrors(ret);
		if(ret != cudaSuccess){
			//printf("Error allocating memory!\n");
			return FAIL;
		}
		
		ret = cudaMemcpy2D(Index_Col_Matrix_d, d_pitch_c, entity.Index_Col_Matrix, sizeof(int)*CodeWord_Size, sizeof(int)*CodeWord_Size, MAX_COLUMN_WEIGHT, cudaMemcpyHostToDevice);
		checkCudaErrors(ret);
		if(ret != cudaSuccess){
			//printf("Error allocating memory!\n");
			return FAIL;
		}
		
		return SUCCESS;
	}
	else{
		/*copy waveform data to GPU terminal*/
		ret = cudaMemcpy(waveform_d, waveform, sizeof(float)*CodeWord_Size, cudaMemcpyHostToDevice);
		checkCudaErrors(ret);
		if(ret != cudaSuccess){
			//printf("Error allocating memory!\n");
			return FAIL;
		}
		
		return SUCCESS;
	}
	
}

bool LDPC_Coding_d::CUDA_Data_callback(LDPC_int *de_info_seq){
	ret = cudaMemcpy(de_info_seq, de_info_seq_d, sizeof(LDPC_int)*CodeWord_Size, cudaMemcpyDeviceToHost);
	checkCudaErrors(ret);
	if(ret == cudaSuccess){
		return SUCCESS;
	}
	else{
		//printf("Memory copy from device to host failed.\n");
		return FAIL;
	}
}



/*The decoded sequence has to copy from GPU to CPU*/
/*optimization:d_pitch_c d_pitch_r CodeWord_Size Info_Size can be substituted by the __constant__ defined values*/
/*can use the global qualfier __restrict__ to avoid pointer alias*/
/*constant value can be accessed CUDA_PARA_SET[0] = CodeWord_Size CUDA_PARA_SET[1] = Info_Size CUDA_PARA_SET[2] = d_pitch_r CUDA_PARA_SET[3] = d_pitch_c*/

__global__ void CUDA_Info_Init(LDPC_Coding_d entity_d, float variance){
	int tid_in_grid = blockDim.x*blockIdx.x+threadIdx.x;
	int count;
	
	/* initialization for r0_d and r1_d matrix */	
	if(tid_in_grid<entity_d.Info_Size){
		for (int j = 0;j < MAX_ROW_WEIGHT;j++){
			entity_d.r0_d[entity_d.d_pitch_r*tid_in_grid/4+j] = -1.0;
			entity_d.r1_d[entity_d.d_pitch_r*tid_in_grid/4+j] = -1.0;
				
		}
	}
	
	if(tid_in_grid<entity_d.CodeWord_Size){
		/* initialization for q0 and q1 matrix */

		for (int i = 0;i < MAX_COLUMN_WEIGHT;i++){
			entity_d.q0_d[entity_d.d_pitch_c*i/4+tid_in_grid] = -1.0;
			entity_d.q1_d[entity_d.d_pitch_c*i/4+tid_in_grid] = -1.0;
			//entity_d.q0_d[entity_d.d_pitch_c*i/4+tid_in_grid] = -1.0;
			//entity_d.q1_d[entity_d.d_pitch_c*i/4+tid_in_grid] = -1.0;
		}
		
		/*Prior Probabilities Calculation*/
		count = 0;
		entity_d.p0_d[tid_in_grid] = 1.0/ (1 + exp(2 * entity_d.waveform_d[tid_in_grid] / variance));
		
		entity_d.p1_d[tid_in_grid] = 1 - entity_d.p0_d[tid_in_grid];
		
		entity_d.p0_init_d[tid_in_grid] = entity_d.p0_d[tid_in_grid];
		entity_d.p1_init_d[tid_in_grid] = entity_d.p1_d[tid_in_grid];
		
		
		while (entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*count/4+tid_in_grid] != -1) {					
			entity_d.q0_d[entity_d.d_pitch_c*count/4+tid_in_grid] = entity_d.p0_d[tid_in_grid];
			entity_d.q1_d[entity_d.d_pitch_c*count/4+tid_in_grid] = entity_d.p1_d[tid_in_grid];
			
			count++;	
		}
		
	}

	
}

/*update the information of check-node(r0, r1)*/

__global__ void LDPC_Decoding_P1(LDPC_Coding_d entity_d){
/********************************************the accessing method is name_of_array[pitch*num_of_row*pitch/4+num_of_col]*******************************************************************************/
	/*can be adapted to two-dimensional thread block*/
	int tid_in_grid = blockIdx.x*blockDim.x*blockDim.y+threadIdx.x+threadIdx.y*blockDim.x;           //the position of each thread in a grid, tid_in_block = threadIdx.x
	int temp_thread = blockIdx.x*blockDim.y+threadIdx.y;
	int outer_count, inner_count, tmp;
	float product;
	/*use shared memory*/
	__shared__ int tmp_row_index_matrix[320];
	
	/*step 1: update the probability of check nodes with probability of bit nodes*/

	if(tid_in_grid < entity_d.Info_Size*blockDim.x){

		tmp_row_index_matrix[threadIdx.y*blockDim.x+threadIdx.x] = entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*temp_thread/4+threadIdx.x];
		__syncthreads();
		
		//while (entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tid_in_grid/4+count] != -1) {						     //update the r0 r1 in the row(th) check formula to count(th) variable node 
		if(tmp_row_index_matrix[threadIdx.y*blockDim.x+threadIdx.x] != -1){
			outer_count = 0;
			product = 1.0;	
				
			__syncthreads();
			//while (entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tid_in_grid/4+outer_count] != -1) {				 //acquire the node message from other variable nodes first				
			while (tmp_row_index_matrix[threadIdx.y*blockDim.x+outer_count] != -1) {	
				inner_count = 0;
				//tmp = entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tid_in_grid/4+outer_count];                     							 //tmp stores the index to search column
				tmp = tmp_row_index_matrix[threadIdx.y*blockDim.x+outer_count]; 
				//while (entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*inner_count/4+tmp] != tid_in_grid) {			 //compare whether two coordinates equal to each other
				while (entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*inner_count/4+tmp] != temp_thread) {	
					inner_count++;
				}
				//if (count != outer_count) {									 			 //exclude the column itself		
				if(threadIdx.x != outer_count) {
				
					product *= (1 - 2*entity_d.q1_d[entity_d.d_pitch_c*inner_count/4+tmp]);			  //a product is responsible for calculating the result of a single thread
	
					__syncthreads();
				}
				
				outer_count++;
			}
			//entity_d.r0_d[entity_d.d_pitch_r*tid_in_grid/4+count] = (1.0 + product) *0.5;
			//entity_d.r1_d[entity_d.d_pitch_r*tid_in_grid/4+count] = (1.0 - product) *0.5;
			tmp = entity_d.d_pitch_r*temp_thread/4+threadIdx.x;
			entity_d.r0_d[tmp] = (1.0 + product) *0.5;
			entity_d.r1_d[tmp] = (1.0 - product) *0.5;

			//count++;
			
		}	
	}
	__syncthreads();
			
} 

/*update the prior credibility(p0, p1) and the information of bit-node(q0, q1)*/
__global__ void LDPC_Decoding_P2(LDPC_Coding_d entity_d){
	int tid_in_grid = blockIdx.x*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x;           //the position of each thread in a grid, tid_in_block = threadIdx.x
	int temp_thread = blockIdx.x*blockDim.y+threadIdx.y;
	int count, inner_count, tmp;
	float product, product2, sum, res1, res2;
	/*use shared memory*/
	__shared__ float p0_sh[DEFAULT_CUDA_THREAD_NUM/4];
	__shared__ float p1_sh[DEFAULT_CUDA_THREAD_NUM/4];
	__shared__ float p0_init_sh[DEFAULT_CUDA_THREAD_NUM/4];
	__shared__ float p1_init_sh[DEFAULT_CUDA_THREAD_NUM/4];
		
	__shared__ int tmp_col_index_matrix[MAX_COLUMN_WEIGHT*DEFAULT_CUDA_THREAD_NUM/4];
	
	/*step 2: update the credibility of bit nodes at iter iteration */
	if(tid_in_grid < entity_d.CodeWord_Size*blockDim.x){
		product = 1.0;
		product2 = 1.0;
		count = 0;
		
		p0_sh[threadIdx.y] = entity_d.p0_d[temp_thread];
		p1_sh[threadIdx.y] = entity_d.p1_d[temp_thread];
		p0_init_sh[threadIdx.y] = entity_d.p0_init_d[temp_thread];
		p1_init_sh[threadIdx.y] = entity_d.p1_init_d[temp_thread];
		
		tmp_col_index_matrix[blockDim.x*threadIdx.y+threadIdx.x] = entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread];		
		__syncthreads();
		
		while (tmp_col_index_matrix[blockDim.x*threadIdx.y+count] != -1) {	
		//while (entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*count/4+temp_thread] != -1) {												//update the p0 p1
			tmp = tmp_col_index_matrix[blockDim.x*threadIdx.y+count];
			//tmp = entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*count/4+temp_thread];
			inner_count = 0;
			//while (entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tmp/4+inner_count] != tid_in_grid) {
			while (entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tmp/4+inner_count] != temp_thread) {
				inner_count++;
			}
			
			product *= entity_d.r0_d[entity_d.d_pitch_r*tmp/4+inner_count];
			product2 *= entity_d.r1_d[entity_d.d_pitch_r*tmp/4+inner_count];
			__syncthreads();
			
			count++;
		}
		
		res1 = __fdividef(p0_init_sh[threadIdx.y] * product, p0_init_sh[threadIdx.y] * product + p1_init_sh[threadIdx.y] * product2);
		res2 = __fdividef(p1_init_sh[threadIdx.y] * product2, p0_init_sh[threadIdx.y] * product + p1_init_sh[threadIdx.y] * product2);
		
		/*set a tolerance 1e-10 to avoid the not-a-number case*/
		//entity_d.p0_d[temp_thread] = (entity_d.p0_d[temp_thread]<1e-10)? 0 : entity_d.p0_init_d[temp_thread] * product / (entity_d.p0_init_d[temp_thread] * product + entity_d.p1_init_d[temp_thread] * product2);
		//entity_d.p1_d[temp_thread] = (entity_d.p0_d[temp_thread]<1e-10)? 1 :entity_d.p1_init_d[temp_thread] * product2 / (entity_d.p0_init_d[temp_thread] * product + entity_d.p1_init_d[temp_thread] * product2);
		
		p0_sh[threadIdx.y] = (p0_sh[threadIdx.y]<1e-10)? 0 : res1;
		p1_sh[threadIdx.y] = (p0_sh[threadIdx.y]<1e-10)? 1 : res2;
				
		__syncthreads();
		
		/*step 3 update the probability of bit nodes with probability of check nodes*/
		if(tmp_col_index_matrix[blockDim.x*threadIdx.y+threadIdx.x] != -1) {
			//tmp = entity_d.Index_Col_Matrix_d[entity_d.d_pitch_c*count/4+tid_in_grid];
			tmp = tmp_col_index_matrix[blockDim.x*threadIdx.y+threadIdx.x];
			inner_count = 0;
			
			//while (entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tmp/4+inner_count] != tid_in_grid) {
			while (entity_d.Index_Row_Matrix_d[entity_d.d_pitch_r*tmp/4+inner_count] != temp_thread) {
				inner_count++;
			}
			
			//entity_d.q0_d[entity_d.d_pitch_c*count/4+tid_in_grid] = (float)entity_d.p0_d[tid_in_grid] / entity_d.r0_d[entity_d.d_pitch_r*tmp/4+inner_count];
			//entity_d.q1_d[entity_d.d_pitch_c*count/4+tid_in_grid] = (float)entity_d.p1_d[tid_in_grid] / entity_d.r1_d[entity_d.d_pitch_r*tmp/4+inner_count];

			entity_d.q0_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread] = (float)p0_sh[threadIdx.y] / entity_d.r0_d[entity_d.d_pitch_r*tmp/4+inner_count];
			entity_d.q1_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread] = (float)p1_sh[threadIdx.y] / entity_d.r1_d[entity_d.d_pitch_r*tmp/4+inner_count];
			
			//sum = entity_d.q0_d[entity_d.d_pitch_c*count/4+tid_in_grid] + entity_d.q1_d[entity_d.d_pitch_c*count/4+tid_in_grid];
			sum = entity_d.q0_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread] + entity_d.q1_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread];
			
			__syncthreads();
			
			//entity_d.q0_d[entity_d.d_pitch_c*count/4+tid_in_grid] /= (float)sum;
			//entity_d.q1_d[entity_d.d_pitch_c*count/4+tid_in_grid] /= (float)sum;
			entity_d.q0_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread] /= (float)sum;
			entity_d.q1_d[entity_d.d_pitch_c*threadIdx.x/4+temp_thread] /= (float)sum;
			
			__syncthreads();
			//count++;
			}
		
		/* step 5 hard decision */ 
		//entity_d.de_info_seq_d[tid_in_grid] = (entity_d.p1_d[tid_in_grid] > 0.5) ? LDPC_1 : LDPC_0;
		entity_d.de_info_seq_d[temp_thread] = (p1_sh[threadIdx.y] > 0.5) ? LDPC_1 : LDPC_0;
		entity_d.p0_d[temp_thread] = p0_sh[threadIdx.y];
		entity_d.p1_d[temp_thread] = p1_sh[threadIdx.y];
		
	}

	__syncthreads();
}




