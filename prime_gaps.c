#include <stdio.h>
#include <stdlib.h>
#include "/usr/lib/x86_64-linux-gnu/mpich/include/mpi.h"

int isPrime(int n) {
		if (n <= 1) return 0;
    if (n == 2 || n == 3) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    
    for (int i=5; i*i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return 0;
    }
    return 1;
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
		
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();

		

		/*
		// check for primes less than N
    // only check numbers of the form (6k + 1) or (6k - 1)
    int count = 2;
    int gap = 2;
    int largest_gap = 2;
    int larger = 5;
    int prev = 5;
    int nums_checked = 2;
    for (int i = 5; i < 100000000; i += 6) {
        if (isPrime(i)) {
            //printf("%d ", i);
            gap = i - prev;
            if (gap > largest_gap) {
                largest_gap = gap;
                larger = i;
            }
            prev = i;
            count++;
        }
        if (isPrime(i+2)) {
            //printf("%d\n", i);
            gap = i+2 - prev;
            if (gap > largest_gap) {
                largest_gap = gap;
                larger = i+2;
            }
            prev = i+2;
            count++;
        }
        nums_checked+=2;
    }
    printf("gap: %d, smaller: %d, larger: %d, count: %d, nums_checked: %d\n", largest_gap, larger-largest_gap, larger, count, nums_checked);
    return 0;
		*/
		
}

