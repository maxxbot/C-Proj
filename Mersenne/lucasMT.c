#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
 
 
 
int lucas_lehmer(unsigned long p)
{
  mpz_t V, mp, t;
  unsigned long k, tlim;
  int res;
 
  if (p == 2) return 1;
  if (!(p&1)) return 0;
 
  mpz_init_set_ui(t, p);
  if (!mpz_probab_prime_p(t, 25)) /* if p is composite, 2^p-1 is not prime */
    { mpz_clear(t); return 0; }
 
  if (p < 23)                     /* trust the PRP test for these values */
    { mpz_clear(t); return (p != 11); }
 
  mpz_init(mp);
  mpz_setbit(mp, p);
  mpz_sub_ui(mp, mp, 1);
 
  /* If p=3 mod 4 and p,2p+1 both prime, then 2p+1 | 2^p-1.  Cheap test. */
  if (p > 3 && p % 4 == 3) {
    mpz_mul_ui(t, t, 2);
    mpz_add_ui(t, t, 1);
    if (mpz_probab_prime_p(t,25) && mpz_divisible_p(mp, t))
      { mpz_clear(mp); mpz_clear(t); return 0; }
  }
 
  /* Do a little trial division first.  Saves quite a bit of time. */
  tlim = p/2;
  if (tlim > (ULONG_MAX/(2*p)))
    tlim = ULONG_MAX/(2*p);
  for (k = 1; k < tlim; k++) {
    unsigned long q = 2*p*k+1;
    /* factor must be 1 or 7 mod 8 and a prime */
    if ( (q%8==1 || q%8==7) &&
         q % 3 && q % 5 && q % 7 &&
         mpz_divisible_ui_p(mp, q) )
      { mpz_clear(mp); mpz_clear(t); return 0; }
  }
 
  mpz_init_set_ui(V, 4);
  for (k = 3; k <= p; k++) {
    mpz_mul(V, V, V);
    mpz_sub_ui(V, V, 2);
    /* mpz_mod(V, V, mp) but more efficiently done given mod 2^p-1 */
    if (mpz_sgn(V) < 0) mpz_add(V, V, mp);
    /* while (n > mp) { n = (n >> p) + (n & mp) } if (n==mp) n=0 */
    /* but in this case we can have at most one loop plus a carry */
    mpz_tdiv_r_2exp(t, V, p);
    mpz_tdiv_q_2exp(V, V, p);
    mpz_add(V, V, t);
    while (mpz_cmp(V, mp) >= 0) mpz_sub(V, V, mp);
  }
  res = !mpz_sgn(V);
  mpz_clear(t); mpz_clear(mp); mpz_clear(V);
  return res;
}


void *mersenne0(void *param0)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne0\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 127; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne1(void *param1)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne1\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 129; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne2(void *param2)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne2\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 131; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne3(void *param3)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne3\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 133; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne4(void *param4)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne4\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 135; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne5(void *param5)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne5\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 137; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne6(void *param6)
{
  clock_t end;
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne6\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 139; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

void *mersenne7(void *param6)
{
  struct timespec start, finish;
  double elapsed;
  unsigned long i;
  double totalPrime;
  printf("Mersenne7\n");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 141; ; i = i + 16) 
  {
    if (lucas_lehmer(i)) 
    {
	  clock_gettime(CLOCK_MONOTONIC, &finish);
	  elapsed = (finish.tv_sec - start.tv_sec);
          elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("M%lu \t%lf\n", i, elapsed);
    }
  }
}

 
int main(int argc, char* argv[]) 
{
  
  pthread_t thread0;
  pthread_t thread1;
  pthread_t thread2;
  pthread_t thread3;
  pthread_t thread4;
  pthread_t thread5;
  pthread_t thread6;
  pthread_t thread7;
  
  
  pthread_create(&thread0, NULL, mersenne0, (void *)thread0);
  pthread_create(&thread1, NULL, mersenne1, (void *)thread1);
  pthread_create(&thread2, NULL, mersenne2, (void *)thread2);
  pthread_create(&thread3, NULL, mersenne3, (void *)thread3);
  pthread_create(&thread4, NULL, mersenne4, (void *)thread4);
  pthread_create(&thread5, NULL, mersenne5, (void *)thread5);
  pthread_create(&thread6, NULL, mersenne6, (void *)thread6);
  pthread_create(&thread7, NULL, mersenne7, (void *)thread7);
  
  while(1)
  {
	  
  }
    
  return 0;
}













































