# TemporalGillespieAlgorithm
C++ code for implementations of the temporal Gillespie algorithm for simulation of epidemic processes on time-varying networks. The algorithm is described in the paper: arXiv:1504.01298 (http://arxiv.org/abs/1504.01298).

The programs need the boost library installed (http://www.boost.org/). 
Compile with the -O2 option for optimal speed, e.g., using g++ the code may be compiled as:

g++ "input name" -o "output name" -O2 -I"path of boost"

All programs simulates independent realizations of the given contagion process on a temporal network given on the form: (t i j) with one triple per line and t,i, and j separated by tabs. ("\t"). 

The programs found here are:
- SIR-Poisson-homogeneous.cpp : Constant-rate (Poissonian) SIR process in a homogeneous population (i.e., same infection/recovery rates for all nodes).
- SIR-Poisson-heterogeneous.cpp : Constant-rate (Poissonian) SIR process in a hoterogeneous population (i.e., infection/recovery rates may differ between nodes).
- SIR-nonMarkovian.cpp : Non-Markovian SIR process with Weibull distributed (https://en.wikipedia.org/wiki/Weibull_distribution) recovery times of individual nodes.
- SIS-Poisson-homogeneous.cpp : Constant-rate (Poissonian) SIS process in a homogeneous population (i.e., same infection/recovery rates for all nodes).
- SIS-Poisson-heterogeneous.cpp : Constant-rate (Poissonian) SIS process in a hoterogeneous population (i.e., infection/recovery rates may differ between nodes).

The Pseudo-Random Number Generator (Mersenne Twister 19937) is initialized with the following seed: 9071982. It can be changed inline. Note however that simulations performed with randomly chosen different seeds are not guaranteed to be (pseudo)independent.

