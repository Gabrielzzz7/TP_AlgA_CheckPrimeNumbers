#include <stdio.h> 
#include <gmp.h> 

int main(void) {
  mpz_t x, y, resultado;


  mpz_init_set_str(x, "2", 10); 
  mpz_init_set_str(y, "6", 10);
  mpz_init(resultado);



  mpz_mul(resultado, x, y); 
  gmp_printf("%Zd \n " "* \n " "%Zd \n " "-------------------- \n " "%Zd \n ", x, y, resultado);


  /* liberar memória usada */
  mpz_clear(x); 
  mpz_clear(y); 
  mpz_clear(resultado);



  return 0;
}