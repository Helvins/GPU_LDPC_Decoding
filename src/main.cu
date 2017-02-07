/*The library in host terminal*/
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
/*The library in device terminal*/
#include <cuda_runtime_api.h>
/*local headers*/
#include "../include/LDPC_Coding.h"
#include "../include/channel.h"

using namespace std;
/*global variables*/
/*
const int cuda_block_num = 256;
const int cuda_thread_num = 256;
*/


int main(int argc, char *argv[]){
	
	cudaEvent_t start, stop;										//create cuda event for recording running time 
	int Info_Size,CodeWord_Size;
	float Code_Rate;
	cudaError_t ret;
	/*set the default number of iteration as 10*/
	int max_iter = (argc>1)? atoi(argv[1]):10;
	int num_set, de_sum = 0;							//num_set stands for the number of decoding sets 
	float variance, error_rate, elapsedTime;
	
	/*initialization of basic parameters*/
#ifdef _LDPC_CODING_H
	Info_Size = DEFAULT_UNCODEWORD_SIZE;
	CodeWord_Size = DEFAULT_CODEWORD_SIZE;
	Code_Rate = DEFAULT_CODE_RATE;
#endif
	//cout.flush();
	LDPC_Coding entity;
	LDPC_Coding_d entity_d;
	Channel *chn = new Channel();
	variance = chn->GetVar(Code_Rate);
	
	/*declaration of information bit sequence and code word sequence*/
	LDPC_int *info_seq = (LDPC_int *)malloc(Info_Size*sizeof(LDPC_int));
	LDPC_int *de_info_seq = (LDPC_int *)malloc(CodeWord_Size*sizeof(LDPC_int));
	LDPC_int *code_word = (LDPC_int *)malloc(CodeWord_Size*sizeof(LDPC_int));
	//float *waveform = (float *)malloc(CodeWord_Size*sizeof(float));
	float *waveform;
	ret = cudaMallocHost((void**)&waveform, CodeWord_Size*sizeof(float));
	
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	/*allocate the memory space in the host and device terminal*/
	if(entity.Memory_Space_Allocation()){
		printf("Host memory allocation success.\n");
	}                   
	if(entity_d.Devide_Memory_Space_Allocation()){
		printf("Device memory allocation success.\n");
	}
	
	/*build check matrix and acquire Index_Col_Matrix*/
	if (entity.Check_Matrix_Construct()) {
		printf("Construct Matrix success\n");
	}
	/*acquire Index_Row_Matrix and free memory space for check matrix*/
	if (entity.Scan_Check_Matrix()) {
		printf("Scan Check matrix success, now we have row and column index matrix\n" );
	}
	/*write row and column index matrix to the output files*/
	entity.WriteData();
	
	/*Memory copy operation from host terminal to device terminal(should include error handling!!!)*/
	//entity_d->CUDA_Memcpy_todev(entity);
	
	entity_d.CUDA_Configuration();
	/*copy memory of index matrix from host terminal to device terminal*/
	entity_d.CUDA_Memcpy_todev(entity, true, waveform);
	
	
	
	for (int i = 0; i < CodeWord_Size; i++) {
		code_word[i] = LDPC_0;
	}
	
	entity.Rand_Seq_Generator(info_seq);
	/*encoding process*/
	if (entity.LDPC_Encoding(info_seq, code_word)) {
		printf("LDPC encoding success\n");
	}
	
	if(argc>2){
		num_set = atoi(argv[2]);
	}
	else{
		printf("Please input the number of sets of code words: ");
		scanf("%d", &num_set);
	}
	
	//cudaEventRecord(start, 0 );	//start recording time
	
	for(int set = 0;set<num_set;set++){
		memset(waveform, 0, sizeof(float)*CodeWord_Size);
		
		/*processed by the channel*/
		chn->BPSK_Modulation(code_word, CodeWord_Size, waveform);
		chn->Channel_Transfer(waveform, CodeWord_Size, Code_Rate, true);
		
		entity_d.CUDA_Memcpy_todev(entity, false, waveform);
	    //calls for GPU kernel function
		{
			CUDA_Info_Init<<<DEFAULT_CUDA_BLOCK_NUM, DEFAULT_CUDA_THREAD_NUM>>>(entity_d, variance);
			/*GPU kernel function*/
			cudaEventRecord(start, 0 );	//start recording time
			for(int iter = 0; iter < max_iter; iter++){
				LDPC_Decoding_P1<<<DEFAULT_CUDA_BLOCK_NUM/2, DEFAULT_CUDA_THREAD_NUM>>>(entity_d);	
				LDPC_Decoding_P2<<<DEFAULT_CUDA_BLOCK_NUM, DEFAULT_CUDA_THREAD_NUM>>>(entity_d);
			}
			cudaEventRecord(stop, 0 );									//stop recording time
			
		}
		
		
		//printf("LDPC decoding success, now try to acquire decoded sequence from GPU terminal\n");
		entity_d.CUDA_Data_callback(de_info_seq);
								
		error_rate = entity.ErrorRate_Check(info_seq, de_info_seq);
		de_sum += error_rate*Info_Size;
		
	}
	//cudaEventRecord(stop, 0 );									//stop recording time
	error_rate = (float)de_sum*1.0 / (Info_Size*num_set);
	printf("The error rate is: %.4f\n", error_rate);
	
	ret = cudaEventElapsedTime( &elapsedTime, start, stop );
	checkCudaErrors(ret);
	if (ret == cudaSuccess){
		printf("The elapsed time for LDPC decoding is %f ms\n", elapsedTime);
	}
	
	/*for test use*/	
	/*
	for(int i = 0;i<10;i++){
		int tmp1 = de_info_seq[i]? 1 : 0;
		int tmp2 = de_info_seq[i+16200]? 1 : 0;
		printf("%d: %d, %d:, %d\n", i, tmp1, i+16200, tmp2);
	}*/
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
		
	return 0;
}
