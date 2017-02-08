# GPU_LDPC_Decoding
Parallel LDPC  decoding in CUDA  environment using NVIDIA vedio card

# Set up the project
1. Type ```git clone -b <branch> https://github.com/Helvins/GPU_LDPC_Decoding ``` to get the selected branch of codes
2. In the root directory, type ```./main```

# Guide
Normally there are 4 optional parameters appended in ./main (number_of_iteration number_of_codeset number_of_gridsize number_of_block)

*number_of_iteration: denotes the number of iterations in LDPC decoding, a large number of iterations will return a relatively low bit error rate. If this input parameter is omitted, the system will adopt 10 iteration as default;

*number_of_codeset: denotes the number of set of undecoded codewords. This will support more precise calculation of error bit rate;

*number_of_gridsize: denotes the size of grid in Stream multiprocessor(SM). If this input parameter is omitted, the system will adopt 256 as default; 

*number_of_block:  denotes the size of block in Stream multiprocessor(SM). If this input parameter is omitted, the system will adopt 256 as default; 
