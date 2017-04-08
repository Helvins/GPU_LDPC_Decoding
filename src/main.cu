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
#include "../include/fileoperation.h"

using namespace std;

int main(int argc, char *argv[]){
	
	cudaEvent_t start, stop;										//create cuda event for recording running time 
	int Info_Size,CodeWord_Size;
	int type_baud_rate;
	float Code_Rate;
	cudaError_t ret;	
	float variance, error_rate, elapsedTime;
	dim3 cfg_para[9];
	
	/*input parameters pretreatment*/
	int num_set, de_sum = 0, max_iter;							//num_set stands for the number of decoding sets 
	/*set the default number of iteration as 10*/
	type_baud_rate = (argc>1)? atoi(argv[1]):0;
	max_iter = (argc>2)? atoi(argv[2]):10;
	num_set = (argc>3)? atoi(argv[3]):1;
	
	/*initialization of basic parameters*/
	switch(type_baud_rate){
		//code word size:64800 code rate: 1/2
		case 0:{
			Info_Size = DEFAULT_UNCODEWORD_SIZE;
			CodeWord_Size = DEFAULT_CODEWORD_SIZE;
			Code_Rate = DEFAULT_CODE_RATE;
		}break;
		//code word size:64800 code rate: 1/2
		case 1:{
			Info_Size = 21600;
			CodeWord_Size = 64800;
			Code_Rate = 1.0/3.0;
		}break;
		//code word size:64800 code rate: 1/2
		case 2:{
			Info_Size = 16200;
			CodeWord_Size = 64800;
			Code_Rate = 1.0/4.0;
		}break;
		//code word size:64800 code rate: 2/3
		case 3:{
			Info_Size = 43200;
			CodeWord_Size = 64800;
			Code_Rate = 2.0/3.0;
		}break;
		//code word size:64800 code rate: 3/4
		case 4:{
			Info_Size = 48600;
			CodeWord_Size = 64800;
			Code_Rate = 3.0/4.0;
		}break;
		//code word size:64800 code rate: 9/10
		case 5:{
			Info_Size = 58320;
			CodeWord_Size = 64800;
			Code_Rate = 9.0/10.0;
		}break;
		//code word size:64800 code rate: 1/4
		case 6:{
			Info_Size = 4050;
			CodeWord_Size = 16200;
			Code_Rate = 1.0/4.0;
		}break;
		//code word size:16400 code rate: 1/2
		case 7:{
			Info_Size = 8100;
			CodeWord_Size = 16200;
			Code_Rate = 1.0/2.0;
		}break;
		//code word size:16200 code rate: 8/9
		case 8:{
			Info_Size = 14400;
			CodeWord_Size = 16200;
			Code_Rate = 8.0/9.0;
		}break;
		default:{
			Info_Size = DEFAULT_UNCODEWORD_SIZE;
			CodeWord_Size = DEFAULT_CODEWORD_SIZE;
			Code_Rate = DEFAULT_CODE_RATE;
		}break;
	}
	//cout.flush();
	LDPC_Coding entity(CodeWord_Size, Code_Rate, type_baud_rate);
	LDPC_Coding_d entity_d(CodeWord_Size, Code_Rate);
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

	
	/*Memory copy operation from host terminal to device terminal(should include error handling!!!)*/
	//entity_d->CUDA_Memcpy_todev(entity);
	
	entity_d.CUDA_Configuration(cfg_para);
	/*copy memory of index matrix from host terminal to device terminal*/
	entity_d.CUDA_Memcpy_todev(entity, true, waveform);
	
	for (int i = 0; i < CodeWord_Size; i++) {
		code_word[i] = LDPC_0;
	}
	
	entity.Rand_Seq_Generator(info_seq);

	/*write unencoded to the output files*/
	//entity.WriteData(info_seq, Info_Size);
	
	/*encoding process*/
	if (entity.LDPC_Encoding(info_seq, code_word)) {
		printf("LDPC encoding success\n");
	}
	
	//cudaEventRecord(start, 0 );	//start recording time
	
	for(int set = 0;set<num_set;set++){
		memset(waveform, 0, sizeof(float)*CodeWord_Size);
		
		/*processed by the channel*/
		chn->BPSK_Modulation(code_word, CodeWord_Size, waveform);
		chn->Channel_Transfer(waveform, CodeWord_Size, Code_Rate);
		
		entity_d.CUDA_Memcpy_todev(entity, false, waveform);
	    //calls for GPU kernel function
		{
			CUDA_Info_Init<<<cfg_para[0], cfg_para[1]>>>(entity_d, variance);
			/*GPU kernel function*/
			
			cudaEventRecord(start, 0 );	//start recording time
			for(int iter = 0; iter < max_iter; iter++){	
				LDPC_Decoding_P1<<<cfg_para[3], cfg_para[4]>>>(entity_d);	
				LDPC_Decoding_P2<<<cfg_para[6], cfg_para[7]>>>(entity_d);
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
	/*write unencoded to the output files*/
	//entity.WriteData(de_info_seq, Info_Size);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
		
	
	
	return 0;
}
