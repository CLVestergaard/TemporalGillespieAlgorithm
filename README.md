# TemporalGillespieAlgorithm
C++ code for implementations of the temporal Gillespie algorithm. The algorithm is described in the paper: arXiv:1504.

The programs need the boost library installed (http://www.boost.org/). 
Compile with the -O2 option for optimal speed, e.g., using g++ the code may be compiled as:

g++ <input name> -o <output name> -O2 -I<path of boost>

All programs simulates independent realizations of the given contagion process on a temporal network given on the form: (t i j) with one triple per line and t,i, and j separated by tabs. ("\t"). 


