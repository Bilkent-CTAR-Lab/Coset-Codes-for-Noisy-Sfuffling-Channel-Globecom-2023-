# Coset-Codes-for-Noisy-Sfuffling-Channel-Globecom-2023-
This library contains MATLAB codes that are employed to generate the numerical results of the paper "A Practical Concatenated Coding Scheme for Noisy
Shuffling Channels with Coset-based Indexing" 
available at: https://ieeexplore.ieee.org/document/10437131
and also at: https://arxiv.org/abs/2410.03578


The main codes are:

BSC_RS_polar_to_run_new_cosets.m

which is used to generate the results over a noisy shuffling channel

and:

BSC_RS_polar_random_sampling_to_run

which is used to generate the results over a noisy shuffling and sampling channels

These main programs call the following functions and libraries:

rel_seq1024.m : 
Contains the  5G standardization unique channel-reliability sequence for a polar code of length 1024. 
The sequences for shorter length codes can be extracted from this sequence. 

build_G.m :
generates the polar transform matrix

bit2gf.m and gf2bit.m: 
Convert data from binary format to GF(2^m) format and from GF(2^m) format to binary format, respectively.

polar_encode.m :
Implements the encoding process for the polar code

polar_dec_BSC :
Implements the successive cancellation decoder for the polar code, when data is received through a binary symmetric channel


 
