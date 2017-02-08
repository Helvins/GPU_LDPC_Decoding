#include "../include/channel.h"
#include "../include/LDPC_Coding.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

Channel::Channel(){																					 //default constructor
	SNR = DEFAULT_SNR;
}

Channel::Channel(float SNR_in){
	SNR = SNR_in;
}

float Channel::GetVar(float Code_Rate) {
	float SNRpbit, N0_uncoded, N0, sigma;
	SNRpbit = pow(10, (SNR / 10));                              //Eb/N0 conversion from dB to decimal
	N0_uncoded = 1.0 / SNRpbit;                                   //default Eb = 1
	N0 = N0_uncoded / Code_Rate;
	sigma = sqrt(N0 / 2);
	return sigma*sigma;											//return the variance of noise

}

float Channel::GaussRand(float Code_Rate) {									//generate a float number that observe the gauss distribution
	static float V1, V2, S;
	static int phase = 0;
	float X;

	if (phase == 0) {
		do {
			float U1 = (float)rand()/RAND_MAX;
			float U2 = (float)rand()/RAND_MAX;

			V1 = 2*U1-1;
			V2 = 2*U2-1;
			S = V1*V1+V2*V2;
		} while (S >= 1 || S == 0);

		X = V1*sqrt(-2*log(S)/S);
	}
	else
		X = V2*sqrt(-2*log(S)/S);

	phase = 1-phase;

	return X*sqrt(GetVar(Code_Rate));
}

void Channel::BPSK_Modulation(LDPC_int *code_word, int length, float *waveform) {
	for (int i = 0;i < length;i++) {
		waveform[i] = (code_word[i] == LDPC_0)? (-1.0):1.0;
	}
	//cout << "BPSK modulationg success." << endl;
}

bool Channel::QPSK_Modulation(LDPC_int *code_word, int length, float *code_word_i, float *code_word_q) {
	int length_out;
	if (length % 2) {
		cout << "The length of the code word has to be even number!" << endl;
		return FAIL;
	}
	else
		length_out = length / 2;

	/* start modulation, the modulated symbol is in I+jQ form */
	for (int i = 0;i < length_out;i++) {
		switch (code_word[2 * i]) {				//output to the i sequence
		case LDPC_0:
			code_word_i[i] = D_QPSK;
			break;
		case LDPC_1:
			code_word_i[i] = -D_QPSK;
			break;
		default:
			return FAIL;
		}
		switch (code_word[2*i+1]){				//output to the q sequence
		case LDPC_0:
			code_word_q[i] = D_QPSK;
			break;
		case LDPC_1:
			code_word_q[i] = -D_QPSK;
			break;
		default:
			return FAIL;
		}
	}
	return SUCCESS;
}

void Channel::Channel_Transfer(float *waveform, int length, float Code_Rate, bool random) {
	for (int i = 0;i < length;i++) {
		if (random){
			waveform[i] += GaussRand(Code_Rate);

		}
		else{
			waveform[i] += 2;
		}
		
	}
	cout << "Channel transfer success." << endl;
}

void Channel::Channel_Transfer(float *code_word_i, float *code_word_q, int length, float Code_Rate) {
	for (int i = 0;i < length;i++) {
		code_word_i[i] += GaussRand(Code_Rate);
		code_word_q[i] += GaussRand(Code_Rate);
	}
	//cout << "Channel transfer success." << endl;
}

void Channel::Parallel_to_Serial(float *code_word_i, float *code_word_q, int length, float *waveform) {
	for (int i = 0;i < length;i++) {
		waveform[2 * i] = code_word_i[i];
		waveform[2 * i+1] = code_word_q[i];
	}
}
