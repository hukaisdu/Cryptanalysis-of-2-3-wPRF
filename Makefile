# Makefile — Cryptanalysis of Two Alternating Moduli Weak PRFs (artifact)
# Reference implementations reproduce Tables 2 and 3; the *_fast variants add
# the solveF3 factorisation, buffer reuse, and multithreading (--threads N).

CXX      ?= g++
CXXFLAGS ?= -std=c++11 -O3 -Wall -Wextra
PTHREAD  ?= -pthread

REF  = wprf_2_3_intersection wprf_2_3_quarter
FAST = wprf_2_3_intersection_fast wprf_2_3_quarter_fast

.PHONY: all ref fast clean
all:  $(REF) $(FAST)
ref:  $(REF)
fast: $(FAST)

# reference (single-threaded, one global RNG)
wprf_2_3_intersection: wprf_2_3_intersection.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

wprf_2_3_quarter: wprf_2_3_quarter.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# optimized + multithreaded  (run: ./wprf_2_3_intersection_fast --seed 1 --threads 8)
wprf_2_3_intersection_fast: wprf_2_3_intersection_fast.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(PTHREAD)

wprf_2_3_quarter_fast: wprf_2_3_quarter_fast.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(PTHREAD)

clean:
	rm -f $(REF) $(FAST)
