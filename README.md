# GPU_LDPC_Decoding
Parallel LDPC  decoding in CUDA  environment using NVIDIA vedio card
notice:the decoding part is combined into a single kernel function with 10 iterations. In the kernel function, parallel threads are called 
twice with 32400 and 64800 in total amount separately.
