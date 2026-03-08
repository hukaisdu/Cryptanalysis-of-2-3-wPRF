# cryptanalysis-of-2-3-wprf

===================================================================

## 1. Overview

**Paper:** **Cryptanalysis of Two Alternating Moduli Weak PRFs. ToSC 2026, Issue 1.**

**Authors:** **Kai Hu, Gregor Leander, Håvard Raddum, Arne Sandrib and Aleksei Udovenko**

The codes of 'wprf_2_3_intersection.cpp', and 'wprf_2_3_quarter.cpp' are used to generate the data of Table 2 and 3, respectively.

===================================================================

## 2. Requirements

**Hardware:** 4GB 

**Software:** g++

**Estimated time:** each instance costs less than 2 hours, most instances are fast

===================================================================

## 3. Environment Installation

Only g++ is needed, so most platforms are suitable to run these codes.
The codes were tested with g++ -std=c++11

===================================================================

## 4. Usage and Experiments

These are codes for the power of 2 case and multiple of 4 case.

Use "g++ wprf_2_3_intersection.cpp -o wprf_2_3_intersection -std=c++11 -O3" to compile the code for experiments to generate the data in Table 2. 

Use "g++ wprf_2_3_quarter.cpp -o wprf_2_3_quarter -std=c++11 -O3" to compile the code for experiments to generate the data in Table 3. 

At the begining of wprf_2_3_intersection.cpp, you can change the parameters:

```
// the length of key and input
const int NN = 40; 
// security
const int lambda = NN / 2.5;
````

NN is the input length, i.e., the "n" in the paper. \lambda is the security level, which can be NN/2 (for aggressive setting) or NN/2.5 (for conservative setting)


In both cpp files, you can change the numbers in Line 14 to change the number of N to in Line 17 to change the factor between 2 and 2.5 to produce the aggressive and conservative cases. 


===================================================================

## 5. License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

===================================================================

## 6. Citation

[Kai Hu, Gregor Leander, Håvard Raddum, Arne Sandrib and Aleksei Udovenko]. "Cryptanalysis of Two Alternating Moduli Weak PRFs". ToSC 2026, Issue 1
Artifact: [hukaisdu/cryptanalysis_2_3_wprf](https://github.com/hukaisdu/Cryptanalysis-of-2-3-wPRF)

===================================================================

## 7. Contact

For questions about this artifact: [[kai.hu@sdu.edu.cn](mailto:kai.hu@sdu.edu.cn)]

===================================================================
