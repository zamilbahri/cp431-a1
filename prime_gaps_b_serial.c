#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#define MAX(A, B) (A > B ? A : B)
#define MIN(A, B) (A < B ? A : B)

void print_mpz(char tag[], mpz_t n, char end[]) {
	printf("%s = ", tag);
	mpz_out_str(stdout, 10, n);
	printf("%s", end);
}

// Calculate the logarithmic integral of x
double li(double x) {
    double result = 0;
    for (double n = 2; n < x; n++) {
        result += 1/log(n);
    }
    return result;
}

// The following are dummy functions that call mpz_functions. These are here so
// 	that I can run a profiler
void get_next_prime(mpz_t dest, mpz_t src) {
	mpz_nextprime(dest, src);
}

int mpz_compare(mpz_t op1, mpz_t op2) {
	return mpz_cmp(op1, op2);
}

void mpz_subtract(mpz_t dest, mpz_t src1, mpz_t src2) {
	mpz_sub(dest, src1, src2);
}

void mpz_assign(mpz_t dest, mpz_t src) {
	mpz_set(dest, src);
}

void find_largest_gap(double num_primes, mpz_t end, mpz_t prime, mpz_t gap, mpz_t prev, mpz_t local_primegap[2], mpz_t first_last_primes[2]) {
	// for (double i = 0; i < num_primes; i++) {
	while (1) {
			
			get_next_prime(prime, prime);

			// break loop if prime > end. Record the last prime encountered
			if (mpz_cmp(prime, end) > 0) {
				mpz_set(first_last_primes[1], prev);
				break;
			}

			// calculate the gap between current and previous prime
			mpz_subtract(gap, prime, prev);

			// if this is the largest gap, then record it and the associated prime
			if (mpz_compare(gap, local_primegap[0]) >= 0) {
				mpz_assign(local_primegap[0], gap);
				mpz_assign(local_primegap[1], prime);
			}

			// update previous prime to current prime
			mpz_assign(prev, prime);
		}
}

int main(int argc, char** argv) {
		int rank=0; int size=1;

		// Initialize mpz_t variables
		mpz_t N;
		int flag_N;

		mpz_init(N);
		mpz_set_ui(N,0);

		flag_N = mpz_set_str(N, argv[1], 10);
		assert (flag_N == 0);

		// timer variables
		clock_t t;

		int flag;
		mpz_t start, end, gap, prime, prev, load;

		// Initialize above mpz_t variables and set initial values
		mpz_init(start);
		mpz_init(end);
		mpz_init(gap);
		mpz_init(prime);
		mpz_init(prev);
		mpz_init(load);

		mpz_set_ui(start,0);
		mpz_set_ui(end,0);
		mpz_set_ui(gap,0);
		mpz_set_ui(prime,2);
		mpz_set_ui(prev,2);
		mpz_set_ui(load,0);
		
		// set the start and end values of this process
		mpz_div_ui(load, N, size);
		mpz_mul_ui(start, load, rank);
		mpz_add(end, start, load);

		double num_primes = li(mpz_get_ui(end));
		//printf("num primes %lf\n", num_primes);

		mpz_t local_primegap[2]; // stores largest gap and the prime associated with it.
		mpz_t first_last_primes[2]; // stores the first and last primes in this process
		unsigned long local_primegap_ui[2];
		unsigned long first_last_primes_ui[2];

		// find the first prime in this process to initialize first_last_primes[]
		get_next_prime(prime, start);
		mpz_set(prev, prime);

		// mpz_init() above arrays
		for (int i = 0; i < 2; ++i) {
			mpz_init(local_primegap[i]);
			mpz_init(first_last_primes[i]);
		}

		mpz_set(local_primegap[0], gap);
		mpz_set(local_primegap[1], prime);

		t = clock();
		find_largest_gap(num_primes, end, prime, gap, prev, local_primegap, first_last_primes);
		t = clock() - t;
		double duration = (double) t/CLOCKS_PER_SEC;

		// get the local largest gap and associated prime in the form of unsigned long
		local_primegap_ui[0] = mpz_get_ui(local_primegap[0]);
		local_primegap_ui[1] = mpz_get_ui(local_primegap[1]);

		// Find the max of time taken by any thread

		// Find the gaps between first primes of a process and last primes of the preceding process
		if (rank == 0) {

			// output
			printf("max gap = %lu, between %lu and %lu\n", local_primegap_ui[0], local_primegap_ui[1]-local_primegap_ui[0], local_primegap_ui[1]);
			printf("global runtime is %lf\n", duration);
		}

		return 0;
}

 