CC = nvcc
CUR_DIR = $(shell pwd)
SRC = $(CUR_DIR)/src
INCLUDE = $(CUR_DIR)/include
BUILD = $(CUR_DIR)/build

CFLAG = -O2
NFLAG = -Wno-deprecated-gpu-targets -m 64 -use_fast_math

OBJ = $(BUILD)/main.o $(BUILD)/LDPC_Coding.o $(BUILD)/channel.o $(BUILD)/fileoperation.o
vpath %.o $(BUILD)

main: main.o LDPC_Coding.o channel.o fileoperation.o
	$(CC) $(CFLAG) $(NFLAG) $(OBJ) -o $@ 

main.o:$(SRC)/main.cu $(INCLUDE)/LDPC_Coding.h $(INCLUDE)/channel.h
	$(CC) $(NFLAG) -c $< -o $(BUILD)/$@

LDPC_Coding.o:$(SRC)/LDPC_Coding.cu $(INCLUDE)/LDPC_Coding.h $(INCLUDE)/fileoperation.h
	$(CC) $(NFLAG) -c $< -o $(BUILD)/$@

channel.o:$(SRC)/channel.cu $(INCLUDE)/channel.h $(INCLUDE)/LDPC_Coding.h
	$(CC) $(NFLAG) -c $< -o $(BUILD)/$@

fileoperation.o:$(SRC)/fileoperation.cu $(INCLUDE)/fileoperation.h
	$(CC) $(NFLAG) -c $< -o $(BUILD)/$@


.PHONY:clean
clean:
	rm $(BUILD)/*.o main
