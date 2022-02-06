#include <stdio.h>
#include <gmp.h>
#include "aux.c"


void numero_aleatorio(mpz_t r, const mpz_t n, gmp_randstate_t rnd) {
    mp_bitcnt_t num_bits = mpz_sizeinbase(n, 2);
    do {
        mpz_urandomb(r, rnd, num_bits);
    } while (!(mpz_cmp_ui(r, 1) >= 0 && mpz_cmp(r, n) <= 0));
}


int talvez_primo(const mpz_t a, const mpz_t n, const mpz_t n1, unsigned int t, const mpz_t q){
  //Caso forem iguais, o teste retorna inconclusivo
  mpz_t fermat_condition, new_a; mpz_inits(fermat_condition, new_a, NULL);
  
  mpz_mod(fermat_condition, a, n);
  mpz_set(new_a, a);
  if(mpz_cmp_ui(fermat_condition, 0) == 0){
    mpz_clears(fermat_condition, new_a, NULL);
    return 1;
  }
  //Caso a for maior que n, basta pegar um número da classe modular de a, t.q 1 <= new_a <= n-1.
  else if(mpz_cmp(a, n) > 0)
  {
    mpz_mod(new_a, a, n);
  }
  mpz_clear(fermat_condition);
  //--------------------------

  mpz_t base, exp, p, two_pow, r; mpz_inits(r, exp, base, p, two_pow, NULL);
  mpz_powm(p, new_a, q, n);
  //Se o número com a menor potência for +- 1, todos os subsequentes serão.
  if(mpz_cmp_ui(p, 1) == 0 || mpz_cmp(p, n1) == 0){
    mpz_clears(base, exp, p, two_pow, r, new_a, NULL);
    return 1;
  }
  mpz_set_ui(two_pow, 2);

  for(unsigned int i = 1; i <= t; i++){
    mpz_pow_ui(r, two_pow, i); mpz_mul(exp, q, r);
    mpz_powm(p, new_a, exp, n);

    //Case ache um -1, todos os outros serão 1's, e não haverá mais alguma contradição para indicar que n é composto.
    if(mpz_cmp(p, n1) == 0){
      mpz_clears(base, exp, p, two_pow, r, new_a, NULL);
      return 1;
    }
    //Diferente de t pois na última iteração se ele for 1 e o algoritmo ainda não tiver parado, deve retornar 1.
    if(mpz_cmp_ui(p, 1) == 0 && i != t){
      mpz_clears(base, exp, p, two_pow, r, new_a, NULL);
      return 0;
    }
  }
  //Se o último número não for 1, falhou no teste de Fermat, e é inconclusivo.
  if(mpz_cmp_ui(p, 1) != 0){
    mpz_clears(base, exp, p, two_pow, r, new_a, NULL);
    return 0;
  }else{
    mpz_clears(base, exp, p, two_pow, r, new_a, NULL);
    return 1;
  }

}



int provavelmente_primo(const mpz_t n, unsigned int iter, gmp_randstate_t rnd){

  mpz_t base, n1, q, aux; mpz_inits(base, n1, q, aux, NULL);
  mpz_sub_ui(n1, n, 1); mpz_set(q, n1);
  unsigned int answer = 0, t = 0;

  //Calcular n-1 = 2^t*q
  while(!mpz_odd_p(q)){
    mpz_tdiv_q_ui(q, q, 2);
    t++;
  }
  mpz_set_ui(aux, 2); mpz_pow_ui(aux, aux, t);
  mpz_tdiv_q(q, n1, aux); mpz_clear(aux);


  for(unsigned int i = 0; i < iter; i++){
    numero_aleatorio(base, n, rnd);
    //Se answer continua igual a i, é porque voltou um 0, e portanto o número é composto.
    answer += talvez_primo(base, n, n1, t, q);

    if(answer == i){
      mpz_clears(base, n1, q, NULL);
      return 0;
    }
  }
  mpz_clears(base, n1, q, NULL);
  return 1;

}

void primo_aleatorio(mpz_t r, unsigned int b, gmp_randstate_t rnd){
  
  unsigned int answer = 0, iter_number = 30;
  while(answer == 0){
    mpz_urandomb(r, rnd, b);
    answer += provavelmente_primo(r, iter_number, rnd);
  }

}
