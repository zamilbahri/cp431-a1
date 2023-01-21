<h1>Prime Gaps</h1>

This program finds the largest gap between two consecutive primes less than some positive integer N. It also outputs two prime numbers in the range [1, N] that realize the largest gap found. Currently, it can run in a reasonable time for N \< 10 billion. 

prime_gaps_a.c only looks at numbers of the form 6k+1 or 6k-1, and the GMP library's mpz_probab_prime_p() function to find primes, and keeps a running track of the gaps to determine the largest one.

prime_gaps_b.c uses the GMP library's mpz_next_prime() function to find primes, and also keeps a running track of the gaps to determine the largest one.

prime_gaps_b_alt.c calculates an estimated number of primes, divides it up across the processors, uses the GMP library's mpz_next_prime() function to find primes, and also keeps a running track of the gaps to determine the largest one.

prime_gaps_serial.c is the same as prime_gaps_b.c, but it doesn't use MPI and only operates on a single processor. I used it for profiling.

<h2>How to run</h2>
It is best to use MacOS or Linux (Ubuntu, for example). The following instructions are for Ubuntu LTS 20.04.

Update Ubuntu:

	sudo apt update && upgrade

Install MPI:

	sudo apt install mpicc

Install GMP:

	sudo apt install libgmp-dev

Compile:

	mpicc -o <output_file>.out <input_file>.c -lmpi -lgmp -lm

Run:

	mpirun -np <number of processors> ./<output_file>.out <N>


