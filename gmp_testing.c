#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char **argv) {
	mpz_t N;
	int flag;

	mpz_init(N);
	mpz_set_ui(N,0);

	flag = mpz_set_str(N, argv[1], 10);
	assert (flag == 0);

	printf("n = ");
	mpz_out_str(stdout, 10, N);
	printf("\n");
	/*
  mpz_t prime;   
  mpz_init(prime);          
  mpz_set_ui(prime, 1);
   
  double i;                                    
  char* num = malloc(sizeof(int)*100000);       
  while(mpz_get_ui(prime) < 100) {
    mpz_nextprime(prime, prime);          
    printf("%lu, \n", mpz_get_ui(prime));   

  }
	*/
}