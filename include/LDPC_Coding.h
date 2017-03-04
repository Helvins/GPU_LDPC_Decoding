#pragma once
#ifndef _LDPC_CODING_H
#define _LDPC_CODING_H

/*macro function definition for error handling in CUDA*/
#define checkCudaErrors(ret) do{if (cudaSuccess != ret) {fprintf(stderr, "Cuda runtime error in line %d, file: %s: %s \n", __LINE__, __FILE__, cudaGetErrorString(ret));break;}}while(0);

typedef bool LDPC_int;
#define LDPC_1 true 
#define LDPC_0 false

#define SUCCESS true
#define FAIL false

#define DEFAULT_CODEWORD_SIZE 64800  														  //size of big blocks for LDPC
#define DEFAULT_CODE_RATE 1.0/2.0
#define DEFAULT_UNCODEWORD_SIZE (unsigned int)(DEFAULT_CODEWORD_SIZE*DEFAULT_CODE_RATE)  	  // namely the size of information bit(32400 in this case)
#define DEFAULT_BLOCK_NUM (unsigned int)(DEFAULT_UNCODEWORD_SIZE/360.0) 					  //(90 in this case)
#define MAX_ROW_WEIGHT 10  																  //the maximum row weight of sparse check matrix
#define MAX_COLUMN_WEIGHT 16																  //the maximum column weight of sparse check matrix, we can get it from the DVB-S2 

/*macro definition for cuda configuration*/
#define DEFAULT_CUDA_BLOCK_NUM 256
#define DEFAULT_CUDA_THREAD_NUM 256


class LDPC_Coding {
public:
	LDPC_Coding();
	LDPC_Coding(int size, float rate);
	~LDPC_Coding();
	bool Memory_Space_Allocation();
	bool Sparse_Cyclic_Matrix_Construct();							  					 	  //Ha matrix
	bool Double_Diagonal_Matrix_Construct();					      						  //Hb matrix
	bool Check_Matrix_Construct();								      						  //combine Ha matrix and Hb matrix together
	bool Calculate_Max_Row_Weight();								  						  //calculate the max row weight in the check matrix
	bool Scan_Check_Matrix();                                         						  //normally return the row index matrix, must be processed after Check_Matrix_Construct()
	void WriteData();
	
	bool LDPC_Encoding(LDPC_int *info_seq, LDPC_int *code_word);      						  //the encoded sequence and uncoded are separately stored in code_word and info_seq
	bool LDPC_Encoding_Check(LDPC_int *code_word);                    						  //Process H*r to verify if the answer equals to vector 0
	//bool LDPC_Decoding(double *waveform, LDPC_int *de_info_seq, int max_iter);	  		  //the decoded sequence and encoded sequence are separately stored in info_seq and code_word
	float ErrorRate_Check(LDPC_int *info_seq, LDPC_int *de_info_seq);						  //return the error bit rate
	/*for test use*/
	void Rand_Seq_Generator(LDPC_int *info_seq);
	
	int *Index_Row_Matrix;											  						  //matrix to store the column position of non-zero elements(Info_Size*MAX_ROW_WEIGHT), the size of each row is the row weight of the sparse matrix 	
	int *Index_Col_Matrix;										      						  //matrix to store the row position of non-zero elements(MAX_COLUMN_WEIGHT*CodeWord_Size), the size of each column is the column weight of the sparse matrix
	
private:
	int CodeWord_Size;												  						  //the size of LDPC codeword, namely N
	float Code_Rate;
	int Info_Size;													  						  //the size of LDPC information sequence, namely K 
	int Blk_Num;													  						  //namely the number of lines in each table
	LDPC_int *Check_Matrix;										  						  	  //check matrix H to supervise whether there are error bit in the codeword
	 
};

/*class for device object*/
class LDPC_Coding_d {
public:
	LDPC_Coding_d();
	LDPC_Coding_d(int size, float rate);
	~LDPC_Coding_d();
	/*allocate memory for p0_d, p1_d, p0_init_d, p1_init_d, q0_d, q1_d, r0_d, r1_d, Index_Col_Matrix_d and Index_Row_Matrix_d*/
	bool Devide_Memory_Space_Allocation();
	/*process dimension and number of thread block and thread*/
    bool CUDA_Configuration(dim3 cfg_para[]);
	/*process data initialization and data memory copy from host terminal to device terminal*/
	
	bool CUDA_Memcpy_todev(LDPC_Coding entity, bool constant, float *waveform);		
	/*copying data of decoded sequence from device terminal to host terminal*/
	bool CUDA_Data_callback(LDPC_int *de_info_seq);
	
	/*data member*/
	float *waveform_d;
	float *p0_d, *p1_d;												  			  //the reliability of x = 1/-1(variable node) in kth iteration
	float *p0_init_d, *p1_init_d;									  			  //the reliability of x = 1/-1(variable node) in the first iteration
	float *q0_d, *q1_d;											 			 	 //the message from variable node to check formula node in kth iteration	
	float *r0_d, *r1_d;											  			 	 //the message from check formula node to variable node in kth iteration
	LDPC_int *de_info_seq_d;											  			  //the decoded sequence stored in GPU device
	
	size_t d_pitch_r;
	size_t d_pitch_c;
	int CodeWord_Size;												  			  //the size of LDPC codeword, namely N
	double Code_Rate;
	int Info_Size;													  			  //the size of LDPC information sequence, namely K 
	int Col_Weight;
	int Row_Weight;
	int *Index_Row_Matrix_d;											  		  //matrix to store the column position of non-zero elements(Info_Size*MAX_ROW_WEIGHT), the size of each row is the row weight of the sparse matrix 	
	int *Index_Col_Matrix_d;										      		  //matrix to store the row position of non-zero elements(MAX_COLUMN_WEIGHT*CodeWord_Size), the size of each column is the column weight of the sparse matrix 
};

/*declaration of kernel function executed in the GPU terminal*/
__global__ void CUDA_Info_Init(LDPC_Coding_d entity_d, float variance);																	  // Initialization for r0 r1 q0 q1
__global__ void LDPC_Decoding_P1(LDPC_Coding_d entity_d);
__global__ void LDPC_Decoding_P2(LDPC_Coding_d entity_d);

#endif // !_LDPC_CODING_H
