#pragma once
#ifndef __CHANNEL_H
#define __CHANNEL_H

#include "LDPC_Coding.h"
#define DEFAULT_SNR 2.0
#define D_QPSK 0.70710678

class Channel{
public:
	Channel();
	Channel(float SNR_in);
	float GetVar(float Code_Rate);
	void BPSK_Modulation(LDPC_int *code_word, int length, float *waveform);
	bool QPSK_Modulation(LDPC_int *code_word, int length, float *code_word_i, float *code_word_q);
	float GaussRand(float Code_Rate);
	void Channel_Transfer(float *waveform, int length, float Code_Rate, bool random);
	void Channel_Transfer(float *code_word_i, float *code_word_q, int length, float Code_Rate);
	void Parallel_to_Serial(float *code_word_i, float *code_word_q, int length, float *waveform);
	
private:
	float SNR;
};



#endif