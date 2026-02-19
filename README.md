These are codes for the power of 2 case and multiple of 4 case.

Use "g++ wprf_2_3_intersection.cpp -o wprf_2_3_intersection -std=c++11 -O3" to compile the code for experiments to generate the data in Table 1. 

Use "g++ wprf_2_3_quarter.cpp -o wprf_2_3_quarter -std=c++11 -O3" to compile the code for experiments to generate the data in Table 2. 

At the begining of wprf_2_3_intersection.cpp, you can change the parameters:

// the length of key and input
const int NN = 40; 
// security
const int lambda = NN / 2.5;

NN is the input length, i.e., the "n" in the paper. \lambda is the security level, which can be NN/2 (for aggressive setting) or NN/2.5 (for conservative setting)


In both cpp files, you can change the numbers in Line 14 to change the number of N to in Line 17 to change the factor between 2 and 2.5 to produce the aggressive and conservative cases. 
# Cryptanalysis-of-2-3-wPRF
