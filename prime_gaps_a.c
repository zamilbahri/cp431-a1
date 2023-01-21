#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <gmp.h>

#define MAX(A, B) (A > B ? A : B)
#define MIN(A, B) (A < B ? A : B)
#define REPS 50

typedef struct ProcessInfo {
	unsigned long first_prime;
	unsigned long last_prime;
} ProcessInfo;

int main(int argc, char** argv) {
    // Initialize the MPI environment
		int rank, size;
		char *ptr;
		unsigned long N = strtoul(argv[1], &ptr, 10);

		double time1, time2, duration, global;

    MPI_Init(&argc, &argv);
		time1 = MPI_Wtime();

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		unsigned long start = 5 + rank * (N/size);
		unsigned long end = MIN(start + (N/size), N);

		mpz_t mpz_prime;
		mpz_init(mpz_prime);

		unsigned long gap = 0; // running gap between primes
		unsigned long prime; // stores the prime associated with the largest gap
		unsigned long prev = 0; // keeps track of previous prime

		unsigned long local_primegap[2]; // stores largest gap and the prime associated with it.
		unsigned long first_last_primes[2]; // stores the first and last primes in this process
		// ProcessInfo pi; // stores the first and last primes in this process

		// Strategy, we only need to check numbers of the form 6k-1 and 6k+1.
		//	Why? Because:
		// 	6k, 6k+2, 6k+4 are divisble by 2
		//	6k is divisble by 3
		//	that leaves 6k+5 (which is same as 6k-1) and 6k+1

		// Increment start until it is of the form 6k-1 or 6k+1
		while ((start % 2 == 0) || (start % 3 == 0)) {
			start++;
		}

		// if start is of the form of 6k+1, then we check separately,
		//	then increment by 4 instead of the usual 6. This makes
		//	the next number of the form 6k-1
		if (start % 6 == 1) {
			mpz_set_ui(mpz_prime, (unsigned long) start);
			if (mpz_probab_prime_p(mpz_prime, REPS) >= 1) {
				first_last_primes[0] = start;
				prime = start;
				prev = start;
			}
			start+=4;
		}

		// Start with number of the form 6k-1, then check
		// 6k+1, then increment by 6
    for (unsigned long i = 0; start + i < end; i+=6) {

				mpz_set_ui(mpz_prime, start+i);
        if (mpz_probab_prime_p(mpz_prime, REPS) >= 1) {
					// set "previous" prime to the first prime encountered in this thread
					if (prev == 0) {
						first_last_primes[0] = start + i;
						prev = start + i;
					}

					if (start + i - prev >= gap) {
						gap = start+i-prev;
						prime = start+i;
					}
					prev = start+i;
				}

				mpz_set_ui(mpz_prime, start + i+2);
        if (start + i+2 < end && mpz_probab_prime_p(mpz_prime, REPS) >= 1) {
					// set "previous" prime to the first prime encountered in this thread
					if (prev == 0) {
						first_last_primes[0] = start + i+2;
						prev = start + i+2;
					}

					if (start + i+2 - prev >= gap) {
						gap = start + i+2 - prev;
						prime = start + i+2;
					}
					prev = start + i+2;
				}
    }
		time2 = MPI_Wtime();
		duration = time2-time1;

		first_last_primes[1] = prev;

		local_primegap[0] = gap;
		local_primegap[1] = prime;
		//printf("rank: %d, start: %d, end: %d, gap: %d prime: %d\n", rank, start, end, gap, prime);

		// Send first prime and last prime info to the root processor (0).
		if (rank != 0) {
			MPI_Send(first_last_primes, 2, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
		}


    unsigned long global_primegap[2];
    MPI_Reduce(local_primegap, global_primegap, 2, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		MPI_Reduce(&duration, &global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);		
		
		// Find if the gaps between the last prime in some thread, and the first prime in the following thread
		//	is larger than the one found using MPI_Reduce()
		if (rank == 0) {

			// initialize an array to store the first and last primes from each processor
			unsigned long all_first_last_primes[size*2];

			// store first and last primes from root processor
			all_first_last_primes[0] = first_last_primes[0];
			all_first_last_primes[1] = first_last_primes[1];

			// get first and last primes from other processors
			MPI_Status status;
			for (int i=1; i < size; ++i) {
				// Can reuse first_last_primes since it'sinfo has already been recorded
				MPI_Recv(first_last_primes, 2, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &status);
				all_first_last_primes[i*2] = first_last_primes[0];
				all_first_last_primes[i*2+1] = first_last_primes;
			}

			// check if previous best prime gap has been beat
			int diff;
			for (int i = 2; i < size*2-1; i+=2) {
				diff = all_first_last_primes[i] - all_first_last_primes[i-1];
				if (diff > global_primegap[0]) {
					global_primegap[0] = diff;
					global_primegap[1] = all_first_last_primes[i];
				}
			}

			printf("Largest gap in primes less than %lu: %lu\n which occured between %lu and %lu\n", N, global_primegap[0], global_primegap[1]-global_primegap[0], global_primegap[1]);
			printf("global runtime is %f\n", global);
		}

    // Finalize the MPI environment.
    MPI_Finalize();
}

 