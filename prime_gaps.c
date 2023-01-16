#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define N 1000
#define MAX(A, B) (A > B ? A : B)
#define MIN(A, B) (A < B ? A : B)

int isPrime(int n) {
		if (n <= 1) return 0;
    if (n <= 3) return 1;
    //if (n % 2 == 0 || n % 3 == 0) return 0;
    
    for (int i=5; i*i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return 0;
    }
    return 1;
}

// float f(float n) {
// 	return (2.0/9)*pow(n, 1.5);
// }

// void findBounds(int n, int size, int rank, float *start, float *end) {
// 	float area = f(n)/size;
	
// 	*start = 0;
// 	for (int i = 0; i <= rank; i++) {
// 		*end = pow((4.5)*(area+f(*start)), 2/3);
// 		if (i != rank) *start = ceil(*end);
// 	}
// 	*start = ceil(*start);
// }

int main(int argc, char** argv) {
    // Initialize the MPI environment
		int rank, size;
		
    MPI_Init(&argc, &argv);

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int start = 5 + rank * (N/size);
		int end = MIN(start + (N/size), N);

		int count = 0; 
		int gap = 0;
		int prime;
		int prev = -1;

		int local_primegap[2];
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
			if (isPrime(start)) {
				prime = start;
				prev = start;
				count++;
			}
			start+=4;
		}

		// Start with number of the form 6k-1, then check
		// 6k+1, then increment by 6
    for (int i = 0; start + i < end; i+=6) {
        if (isPrime(start + i)) {
					// set "previous" prime to the first prime encountered in this thread
					if (prev < 0) prev = start + i;

					if (start + i - prev >= gap) {
						gap = start+i-prev;
						prime = start+i;
					}
					prev = start+i;
					count++;
				}

        if (start + i+2 < end && isPrime(start + i+2)) {
					// set "previous" prime to the first prime encountered in this thread
					if (prev < 0) prev = start + i+2;

					if (start + i+2 - prev >= gap) {
						gap = start + i+2 - prev;
						prime = start + i+2;
					}
					prev = start + i+2;
					count++;
				}
    }

		local_primegap[0] = gap;
		local_primegap[1] = prime;
		//printf("rank: %d, start: %d, end: %d, gap: %d prime: %d\n", rank, start, end, gap, prime);

    int global_primegap[2];
    MPI_Reduce(local_primegap, global_primegap, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);


    if (rank == 0) {
			printf("Largest gap in primes less than %d: %d\n which occured between %d and %d\n", N, global_primegap[0], global_primegap[1]-global_primegap[0], global_primegap[1]);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}

