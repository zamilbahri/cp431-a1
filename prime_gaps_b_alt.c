#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gmp.h>
#include <math.h>
#include <assert.h>

#define MAX(A, B) (A > B ? A : B)
#define MIN(A, B) (A < B ? A : B)
#define ROOT 0

void print_mpz(char tag[], mpz_t n, char end[]) {
	printf("%s = ", tag);
	mpz_out_str(stdout, 10, n);
	printf("%s", end);
}

double li(double x) {
    double result = 0;
    for (double n = 2; n < x; n++) {
        result += 1/log(n);
    }
    return result;
}

double nth_prime(double n) {
	if (n == 0) return 1;
	return n * log(n);
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
		int rank, size;

		// Initialize mpz_t variables
		mpz_t N;
		int flag_N;
		MPI_Status status;

		mpz_init(N);
		mpz_set_ui(N,0);

		flag_N = mpz_set_str(N, argv[1], 10);
		assert (flag_N == 0);

		// timer variables
		double time1, time2, duration, global_duration;

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

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// double *num_primes_less_than = malloc(size * sizeof(double));
		// double *tmp = malloc(size * sizeof(double));
		// if(num_primes_less_than==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
		
		double num_primes; // how many primes I want to loop over in this process
		if (rank == ROOT) {
			num_primes = li(mpz_get_ui(N))/size;
			// printf("num_primes = %lf\n", num_primes);
			for (int i = 1; i < size; ++i) {
				MPI_Send(&num_primes, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		} else {
			MPI_Recv(&num_primes, 1, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD, &status);
		}

		// set the start and end values of this process
		mpz_div_ui(load, N, size);
		mpz_mul_ui(start, load, rank);
		mpz_add(end, start, load);

		// num_primes_less_than[rank] = li(mpz_get_ui(end));
		// MPI_Allgather(&num_primes_less_than[rank], 1, MPI_DOUBLE, tmp, 1, MPI_DOUBLE, MPI_COMM_WORLD);

		// double num_primes = tmp[rank] - (rank == 0 ? 0 : tmp[rank-1]);


		// for (int i = 0; i < size; ++i) {
		// 	if (rank == i) {
		// 		printf("num_primes less than %lu: %lf\n", mpz_get_ui(end), tmp[i]);
		// 		printf("rank: %d num_primes: %lf\n", i, num_primes);
		// 	}
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// }

		// free(num_primes_less_than);
		

		mpz_t local_primegap[2]; // stores largest gap and the prime associated with it.
		mpz_t first_last_primes[2]; // stores the first and last primes in this process
		unsigned long local_primegap_ui[2];
		unsigned long first_last_primes_ui[2];

		// find the first prime in this process to initialize first_last_primes[]
		// unsigned long estimated_startprime = (rank == 0 ? 0 : floor(tmp[rank-1]) * log(floor(tmp[rank-1])));
		// printf("rank:  %d, estimated start prime: %lu\n", rank, estimated_startprime);
		// mpz_set_ui(start, estimated_startprime);

		// unsigned long estimated_start_prime = floor(nth_prime(num_primes * rank));
		// mpz_set_ui(start, estimated_start_prime);

		mpz_nextprime(prime, start);
		mpz_set(prev, prime);

		// free(tmp);

		// mpz_init() above arrays
		for (int i = 0; i < 2; ++i) {
			mpz_init(local_primegap[i]);
			mpz_init(first_last_primes[i]);
		}

		mpz_set(local_primegap[0], gap);
		mpz_set(local_primegap[1], prime);

		mpz_set(first_last_primes[0], prime);

		// start timer
		time1 = MPI_Wtime();

		// for (double i = 0; i < num_primes; ++i) {

		for (double i = 0; i < num_primes; ++i) {
			
			mpz_nextprime(prime, prime);

			// break loop if prime > end. Record the last prime encountered
			// if (mpz_cmp(prime, end) > 0) {
			// 	mpz_set(first_last_primes[1], prev);
			// 	break;
			// }

			// calculate the gap between current and previous prime
			mpz_sub(gap, prime, prev);

			// if this is the largest gap, then record it and the associated prime
			if (mpz_cmp(gap, local_primegap[0]) >= 0) {
				mpz_set(local_primegap[0], gap);
				mpz_set(local_primegap[1], prime);
			}

			// update previous prime to currrent prime
			mpz_set(prev, prime);
		}

		mpz_t leftover; mpz_init(leftover);
		mpz_sub(leftover, end, prev);
		// printf("rank=%d, start=%lu, prev=%lu, end=%lu, leftover=%lu\n", rank, mpz_get_ui(start), mpz_get_ui(prev), mpz_get_ui(end), mpz_get_ui(leftover));

		MPI_Barrier(MPI_COMM_WORLD);

		if (mpz_cmp(end, prev) > 0) {
			// printf("rank %d needs to continue\n", rank);
			while (1) {
				mpz_nextprime(prime, prime);

				// break loop if prime > end. Record the last prime encountered
				if (mpz_cmp(prime, end) > 0) {
					mpz_set(first_last_primes[1], prev);
					break;
				}

				// calculate the gap between current and previous prime
				mpz_sub(gap, prime, prev);

				// if this is the largest gap, then record it and the associated prime
				if (mpz_cmp(gap, local_primegap[0]) >= 0) {
					mpz_set(local_primegap[0], gap);
					mpz_set(local_primegap[1], prime);
				}

				// update previous prime to currrent prime
				mpz_set(prev, prime);
			}
		}


		mpz_set(first_last_primes[1], prev);

		// end timer
		time2 = MPI_Wtime();
		duration = time2-time1;
		// printf("rank = %d, local duration = %lf\n", rank, duration);

		// send first and last primes to root thread
		first_last_primes_ui[0] = mpz_get_ui(first_last_primes[0]);
		first_last_primes_ui[1] = mpz_get_ui(first_last_primes[1]);
		// for (int i = 0; i < size; ++i) {
		// 	if (rank == i){
		// printf("rank: %d, first: %lu, last1: %lu, leftover: %lu\n", rank, first_last_primes_ui[0], first_last_primes_ui[1], mpz_get_ui(leftover));
		// 	}
		// }
		if (rank != 0) {
			MPI_Send(first_last_primes_ui, 2, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
		}

		// get the local largest gap and associated prime in the form of unsigned long
		local_primegap_ui[0] = mpz_get_ui(local_primegap[0]);
		local_primegap_ui[1] = mpz_get_ui(local_primegap[1]);
		
		// Find the max of all the local largest gaps
		unsigned long global_primegap[2];
		MPI_Reduce(local_primegap_ui, global_primegap, 2, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		// Find the max of time taken by any thread
		MPI_Reduce(&duration, &global_duration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		// Find the gaps between first primes of a process and last primes of the preceding process
		if (rank == 0) {
			
			// array to store the first and last primes
			unsigned long all_first_last_primes[size*2];

			// store the first and last primes from root process
			all_first_last_primes[0] = first_last_primes_ui[0];
			all_first_last_primes[1] = first_last_primes_ui[1];

			// receive first and last primes from other processes and store them
			for (int i = 1; i < size; ++i) {
				MPI_Recv(first_last_primes_ui, 2, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &status);
				all_first_last_primes[i*2] = first_last_primes_ui[0];
				all_first_last_primes[i*2+1] = first_last_primes_ui[1];
			}

			// find the difference between consecutive elements (ignore first and last elements).
			//	If diff > largest gap, then recognize it
			unsigned long diff;
			for (int i = 2; i < size*2-1; i+=2) {
				if (all_first_last_primes[i] < all_first_last_primes[i-1]) continue;
				diff =  all_first_last_primes[i] - all_first_last_primes[i-1];
				if (diff > global_primegap[0]) {
					global_primegap[0] = diff;
					global_primegap[1] = all_first_last_primes[i];
				}
			}

			// output
			printf("max gap = %lu, between %lu and %lu\n", global_primegap[0], global_primegap[1]-global_primegap[0], global_primegap[1]);
			printf("global runtime is %f\n", global_duration);
		}

		// End MPI Program
    MPI_Finalize();

		return 0;
}

 