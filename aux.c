#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"

#define STANDARD_SIZE  1024

void mdc_estendido(mpz_t g, mpz_t x, mpz_t y, const mpz_t a, const mpz_t b){

  mpz_t first_num, second_num, rest; mpz_init(first_num); mpz_init(second_num); mpz_init(rest);


  //Place the greatest number always on the 'first_number' variable.
  if(mpz_cmp(a,b) >= 0 ){
   mpz_set(first_num,a);mpz_set(second_num,b);
  }else{
   mpz_set(first_num, b);mpz_set(second_num, a);
  }

  mpz_t quocient[STANDARD_SIZE];

  long long unsigned int loop_number = 0;

  while(mpz_cmp_ui(second_num,0) != 0){
    mpz_init(quocient[loop_number]); mpz_tdiv_q(quocient[loop_number], first_num, second_num);
    mpz_mod(rest, first_num, second_num);
    mpz_set(first_num, second_num); mpz_set(second_num, rest);
    loop_number++; 
  }
  mpz_set(g, first_num); mpz_set_ui(x, 1);mpz_set_ui(y, 0);
  mpz_t new_x, new_y, qz; mpz_init(qz); mpz_init(new_x); mpz_init(new_y);

  //Calculate the X and  Y.
  for(long long unsigned int i = loop_number; i > 0; i--){
    mpz_set(new_x, y); mpz_mul(qz, quocient[i-1], y); mpz_sub(new_y, x, qz);
    mpz_set(x, new_x);mpz_set(y, new_y);mpz_clear(quocient[i-1]);
  }

  if(mpz_cmp(a,b) < 0){
    mpz_set(new_x, y); mpz_set(y, x); mpz_set(x, new_x);
  }
  mpz_clears(first_num, second_num, new_x, new_y, qz, rest, NULL);
}

int inverso_modular(mpz_t r, const mpz_t a, const mpz_t n){

  mpz_t x, y, g; mpz_inits(x, y, g, NULL);
  mdc_estendido(g, x, y, a, n);
  if (mpz_cmp_ui(g, 1) != 0){
    mpz_clears(x, y, g, NULL);

    return 0;
  }else{
    mpz_set(r, g);
    mpz_clears(x, y, g, NULL);
    return 1;
  }
  
}

void exp_binaria(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t n){

  mpz_t base, expoente, temp2; mpz_inits(base, expoente,temp2, NULL);

  mpz_set(expoente, e); mpz_set(base, b);

  mpz_mod(base, base, n);
  if(mpz_cmp_ui(base, 0) == 0){
    mpz_set_ui(r, 0); 
    return;
  }
  mpz_set_ui(r, 1); 
    
  while(mpz_cmp_ui(expoente, 0) > 0){
    
    mpz_mod_ui(temp2, expoente, 2);
    if (mpz_cmp_ui(temp2, 1) == 0){
      mpz_mul(r, r, base);
      mpz_mod(r, r, n);
    }
    mpz_tdiv_q_ui(expoente, expoente, 2);

    mpz_mul(base, base, base); mpz_mod(base, base, n);

  }
  mpz_clears(base, expoente, temp2, NULL);

}

// int main(int argc, char* argv[]){


//   mpz_t a, b, g,x,y;
//   mpz_init(a); mpz_init(g);mpz_init(b);mpz_init(x);mpz_init(y);
//   mpz_set_ui(a, 1100);mpz_set_ui(b, 1282); 
//   mdc_estendido(g, x, y, a, b);
//   gmp_printf("MDC de %Zd e %Zd é: %Zd\nAssim:\n\
//     %Zd*%Zd + %Zd*%Zd = %Zd\n", b, a, g, a, x, b, y, g);

// }
