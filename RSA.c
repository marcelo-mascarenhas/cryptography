#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>

#include "miller_rabin.c"



void gera_chaves(mpz_t n, mpz_t e, mpz_t d, gmp_randstate_t rnd){
  long long unsigned int bit_number = 2048, e_value = 65537;
  mpz_t p, q, fi_n, aux1, aux2, aux3; mpz_inits(p, q, fi_n, aux1, aux2, aux3, NULL);
  
  //Gera p e q.
  primo_aleatorio(p, bit_number, rnd); primo_aleatorio(q, bit_number, rnd);
  //Gera n.
  mpz_mul(n, p, q);   
  //Gera fi_n
  mpz_sub(aux1, p, 1); mpz_sub(aux2, p, 1);  
  mpz_mul(fi_n, aux1, aux2);  mpz_clears(aux1, aux2, NULL);

  //Calcula um valor de 'e' que seja relativamente primo em comparação à fi(n)
  while (1){
    mpz_set_ui(e, e_value);
    mpz_gcd(aux3, fi_n, e);
    if(mpz_cmp_ui(aux3, 1) == 0){
      break;
    }
    e_value+=2;
  }
  inverso_modular(d, e, fi_n);
  mpz_clears(aux3, p, q, fi_n, NULL);

}

void codifica(mpz_t r, const char *str){

  mpz_t aux1, aux2, aux3, base; mpz_inits(aux1, aux2, aux3, base, NULL);
  
  mpz_set_ui(r, 0); mpz_set_ui(base, 256);
  for(unsigned int i = 0; i < strlen(str); i++ ){
    mpz_pow_ui(aux1, base, i); mpz_set_ui(aux2, str[i]);
    mpz_mul(aux3, aux1, aux2); mpz_add(r, r, aux3);
  }
  mpz_clears(aux1, aux2, aux3, base, NULL);
}

char *decodifica(const mpz_t n){
  char* str;
  str = (char*)calloc(sizeof(char), 501);
  mpz_t new_n, base, aux1,aux2, aux3, aux4; mpz_inits(new_n, base, aux1, aux2, aux3, aux4, NULL);
  mpz_set(new_n, n); mpz_set_ui(base, 256); unsigned int i = 1;
  while (mpz_cmp_ui(new_n, 100) >= 0){
    mpz_pow_ui(aux1, base, i); mpz_mod(aux2, n, aux1);
    if(i == 1){
      str[i-1] = mpz_get_ui(aux2);
      mpz_sub(new_n, new_n, aux2); i++;
      continue;
    }
    mpz_pow_ui(aux3, base, i-1); mpz_tdiv_q(aux4, aux2, aux3);
    str[i-1] = mpz_get_ui(aux4);
    mpz_sub(new_n, new_n, aux2);
    i++;
  }
  return str;
}

void criptografa(mpz_t C, const mpz_t M, const mpz_t n, const mpz_t e){
  mpz_powm(C, M, e ,n);
}


void descriptografa(mpz_t M, const mpz_t C, const mpz_t n, const mpz_t d){
  mpz_powm(M, C, d, n );
}
