#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>

#include "miller_rabin.c"


void gera_primo_seguro(mpz_t p, unsigned int b, gmp_randstate_t rnd){

  unsigned int loop_num = 2, greater_loop_num = 5, trueness_p = 0, trueness_q = 0;
  mpz_t q; mpz_init(q);
  while (1){
    //Seleciona um número aleatório e calcula p = 12K+11
    mpz_urandomb(p, rnd, b-7);
    mpz_mul_ui(p, p, 12); mpz_add_ui(p, p, 11);
    //Faz a verificação em cadeia.
    trueness_p = provavelmente_primo(p, loop_num, rnd);
    
    if(trueness_p == 1){
      mpz_sub_ui(q, p, 1); mpz_div_ui(q, q, 2);

      trueness_q = provavelmente_primo(q, loop_num, rnd);
      
      if(trueness_q == 1){
        trueness_p = provavelmente_primo(p, greater_loop_num, rnd);
        
        if(trueness_p == 1){
          trueness_q = provavelmente_primo(q, greater_loop_num, rnd);
          if(trueness_q == 1){
            mpz_clear(q);
            break;
          }
        }
      } 
    }
  }
}


void encontra_gerador(mpz_t g, const mpz_t p, gmp_randstate_t rnd){
  mpz_t q, new_p, aux2, aux3, rnd_i; mpz_inits(rnd_i, q, new_p, aux2, aux3, NULL);
  //Calcular q.
  mpz_set(q, p); mpz_sub_ui(q, q, 1); mpz_div_ui(q, q, 2);
  mpz_set(new_p, p); mpz_sub_ui(new_p, new_p, 1);

  mpz_set_ui(g, 1);
  while(1){
    //Checando as ordens do g aleatório. Se g^2 mod p != 1 (ord 1/2), se g^2 mod p != 1(ord q).
    numero_aleatorio(rnd_i, new_p, rnd);
    mpz_powm(aux2, rnd_i, q, p); mpz_powm_ui(aux3, rnd_i, 2, p);
    if(mpz_cmp_ui(aux2, 1) != 0 && mpz_cmp_ui(aux3, 1) !=0){
      mpz_set(g, rnd_i);
      break;
    }
  }

}

void gera_chave_elgamal(mpz_t p, mpz_t g, mpz_t A, gmp_randstate_t rnd){
  unsigned long int b = 512;
  gera_primo_seguro(p, b, rnd); encontra_gerador(g, p, rnd);
  mpz_t a, new_p; mpz_inits(a, new_p, NULL); mpz_sub_ui(new_p, p, 1);
  numero_aleatorio(a, new_p, rnd);
  mpz_powm(A, g, a, p);
}

void criptografa(mpz_t B, mpz_t C, const mpz_t M, const mpz_t p, const mpz_t g, const mpz_t A, gmp_randstate_t rnd){

  //Gera o g^b aleatório.
  mpz_t b, new_p, aux; mpz_inits(aux, b, new_p, NULL); mpz_sub_ui(new_p, p, 1);
  numero_aleatorio(b, new_p, rnd);
  mpz_powm(B, g, b, p);
  //Calcula A^b.
  mpz_powm(aux, A, b, p);
  mpz_mul(C, aux, M);

}

void descriptografa(mpz_t M, const mpz_t B, const mpz_t C, const mpz_t p, const mpz_t g, const mpz_t a){
  mpz_t aux, inverse_mod; mpz_inits(aux, inverse_mod, NULL);
  mpz_powm(aux, B, a, p);
  mpz_invert(inverse_mod, aux, p);
  mpz_mul(M, inverse_mod, C);
  mpz_mod(M, M, p);
}

int main(){
  mpz_t g, p; mpz_inits(g,p, NULL); mpz_set_ui(p, 357655694505827); mpz_set_ui(g, 1);
  long long int abc = 59843209543;
  gmp_randstate_t rnd; gmp_randinit_mt(rnd); gmp_randseed_ui(rnd, abc);
  encontra_gerador(g, p, rnd);
}
