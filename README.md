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

**Estimated time:**, when NN is small, the instances are fast. When NN is larger such as 64, the codes may run for hours. 

===================================================================

## 3. Environment Installation

Only g++ is needed, so most platforms are suitable to run these codes.
The codes were tested with g++ -std=c++11

===================================================================

## 4. Usage and Experiments

These are codes for the power of 2 case and multiple of 4 case.

Use `g++ wprf_2_3_intersection.cpp -o wprf_2_3_intersection -std=c++11 -O3` to compile the code for experiments to generate the data in Table 2. 

Use `g++ wprf_2_3_quarter.cpp -o wprf_2_3_quarter -std=c++11 -O3` to compile the code for experiments to generate the data in Table 3. 

Alternatively, a `Makefile` is provided:

```
make        # build all four binaries
make ref    # only the reference binaries (Tables 2 and 3)
make fast   # only the optimized multithreaded variants (see Section 4.1)
make clean
```

At the begining of wprf_2_3_intersection.cpp, you can change the parameters:

```
// the length of key and input
const int NN = 32;
// security
const int lambda = NN / 2;
```

NN is the input length, i.e., the "n" in the paper. \lambda is the security level, which can be NN/2 (for aggressive setting) or NN/2.5 (for conservative setting)

To change the parameters, edit these lines directly:

- `wprf_2_3_intersection.cpp`: `NN` on **Line 18**, `lambda` on **Line 21**.
- `wprf_2_3_quarter.cpp`: `NN` on **Line 16**, `lambda` on **Line 19**.

Change `NN` to set the input length `n`, and change the factor in `lambda = NN / 2` (aggressive) or `lambda = NN / 2.5` (conservative) to switch between the two settings.

### Mapping the program output to the paper tables

Each run prints the log2 of the average (over 100 trials) of several counters. They map to the columns of Table 2 / Table 3 as follows:

| Program output        | Column in Table 2 / Table 3 |
|-----------------------|-----------------------------|
| `Sampling1`, `Sampling2` | "Sampling x" (the two are independent sampling counters; the table reports their average `(Sampling1 + Sampling2) / 2`) |
| `W1_before`           | `#w` before the {0,1}-filter |
| `W1`, `W2`            | `#w` after the {0,1}-filter (the two agree up to noise; the table reports a single value) |
| `Exhaustive`          | `#keys from w` (candidates left for the final exhaustive check) |
| `Success`             | number of successful runs out of 100 (not a table column) |

### Running with a fixed seed (reproducibility)

Each program draws a fresh nondeterministic seed from `std::random_device` on every run, so the exact per-run figures vary slightly. To reproduce a specific run, pass `--seed N`:

```
./wprf_2_3_intersection --seed 1
./wprf_2_3_quarter --seed 1
```

The seed actually used is printed at the top of the output. With the same seed, the same parameters, and the same binary, a run is reproducible.

`wprf_2_3_intersection.cpp` has a `verbose` flag (Line 15), set to `false` by default, so only the parameter line, a progress bar, and the final summary are printed; set it to `true` for the detailed per-trial trace.

### Success rate

The `Success` line counts, out of the 100 trials, how many recovered the correct key (in `wprf_2_3_intersection`) / the correct `w'` (in `wprf_2_3_quarter`). The quantities reported in Tables 2 and 3 (the Sampling, `#w`, and `#keys` columns) are averages over all 100 trials and do not depend on this count.

A success rate below 100% is expected for some parameter sets: it comes from the randomized sampling occasionally failing to produce, within a single run, an input pair of the required differential rank. It is not a fixed theoretical figure. Key recovery itself is deterministic once such a pair is found, because the surviving key candidates are verified against known input/output pairs.

### 4.1 Optimized, multithreaded variants

`wprf_2_3_intersection_fast.cpp` and `wprf_2_3_quarter_fast.cpp` are faster drop-in variants that reproduce the same Tables 2 and 3. They change nothing in the attack itself; they add:

- **An algorithmic optimization.** In the `#w(before filter)` loop the F3 coefficient matrix is constant across the inner guess loop, so it is factored once and only back-substituted for each guess (plus scratch-buffer reuse). This produces exactly the same solution sets as the reference `solveF3`.
- **Multithreading.** The 100 trials are independent and are distributed over threads. Pass `--threads N` (default `1`):

```
make fast
./wprf_2_3_quarter_fast --seed 1 --threads 8
```

Each trial is seeded from `(seed, trial_index)`, so the output is **independent of the thread count**: `--threads 1` and `--threads N` give bit-for-bit identical results for the same `--seed`. Because this per-trial seeding differs from the reference programs' single global generator, the exact per-run numbers differ slightly from the reference binaries, but the averages match the same theory and experiments in Tables 2 and 3.

The reference `wprf_2_3_intersection.cpp` / `wprf_2_3_quarter.cpp` remain the primary artifacts for the tables; the `_fast` variants are provided only to make the larger parameter sets (e.g. `n = 64`, or the conservative `n = 44` case) practical to run.


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
