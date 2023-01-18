<h1>Prime Gaps</h1>

This program finds the largest gap between two xonsecutive primes less than N. It also outputs two prime numbers in the range [1, N] that realize the largest gap found. Currently, it can run in a reasonable time for N \< 1 billion. 

prime_gaps_a.c uses the sieve of Eratosthenes, and the GMP library's mpz_probab_prime_p() function to find primes, and keeps a running track of the gaps to determine the largest one.

prime_gaps_b.c uses the GMP library's mpz_next_prime() function to find primes, and also keeps a running track of the gaps to determine the largest one.

<h2>Known Issues</h2>

* Takes roughly 6 seconds on my laptop for N = 100 million for both methods, although prime_gaps_a.c is slightly faster. Need to scale better for 1 billion and maybe 1 trillion.

