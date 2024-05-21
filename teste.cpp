#include <iostream>
#include <gmp.h>

// a principio, o codigo funciona para grande parte dos problemas pedidos. A ultima parte ficou problematica (3 questao)

bool millerRabinTest(mpz_t n, mpz_t a) {
  mpz_t x, n_minus_1, d, r;
  mpz_inits(x, n_minus_1, d, r, NULL);

  // Calcula n - 1
  mpz_sub_ui(n_minus_1, n, 1);

  // Inicializa d como n - 1
  mpz_set(d, n_minus_1);

  // Calcula o valor de d e r tal que d * 2^r = n - 1
  mpz_set_ui(r, 0);
  while (mpz_even_p(d)) {
    mpz_divexact_ui(d, d, 2);
    mpz_add_ui(r, r, 1);
  }

  // Calcula x = a^d % n
  mpz_powm(x, a, d, n);

  // Se x == 1 ou x == n - 1, n é provavelmente primo
  if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) {
    mpz_clears(x, n_minus_1, d, r, NULL);
    return true;
  }

  // Realiza r - 1 iterações
  mpz_t i;
  mpz_init_set_ui(i, 0);
  for (i;; mpz_add_ui(i, i, 1)) {
    if (mpz_cmp(i, r) >= 0) break;

    // Calcula x = x^2 % n
    mpz_powm_ui(x, x, 2, n);

    // Se x == n - 1, n é provavelmente primo
    if (mpz_cmp(x, n_minus_1) == 0) {
      mpz_clears(x, n_minus_1, d, r, NULL);
      mpz_clear(i);
      return true;
    }

    // Se x == 1, n é composto
    if (mpz_cmp_ui(x, 1) == 0) {
      mpz_clears(x, n_minus_1, d, r, NULL);
      mpz_clear(i);
      return false;
    }
  }

  mpz_clears(x, n_minus_1, d, r, NULL);
  return false; // n é composto
}

// Teste de primalidade de Miller-Rabin
bool isPrimeMRtest(mpz_t n, int& MR_count) {
  mpz_t primes[40];

  // Inicializa os números primos
  const char* primes_str[40] = {
      "2", "3", "5", "7", "11", "13", "17", "19", "23", "29",
      "31", "37", "41", "43", "47", "53", "59", "61", "67", "71",
      "73", "79", "83", "89", "97", "101", "103", "107", "109", "113",
      "127", "131", "137", "139", "149", "151", "157", "163", "167", "173"
  };

  for (int i = 0; i < 40; ++i) {
    mpz_init_set_str(primes[i], primes_str[i], 10);
  }

  bool isPrime = true;

  for (int i = 0; i < 40; i++) {
    MR_count++;
    if (!millerRabinTest(n, primes[i])) { // Usa os 40 primeiros primos
      isPrime = false;
      break;
    }
  }

  for (int i = 0; i < 40; ++i) {
    mpz_clear(primes[i]);
  }

  if (isPrime) {
    return true; // n é primo
  }

  return false; // n é composto
}



int main() {
  mpz_t N, a, candidate, prime, generator;
  int MR_count = 0;

  mpz_init(N);
  mpz_init(a);
  mpz_init(candidate);
  mpz_init(prime);
  mpz_init(generator);


  // Leitura
  std::string input_N, input_a;
  std::cout << "Enter a large number N: ";
  // printf("ssssssssssssss\n");
  std::cin >> input_N;
  // printf("adsdaadads\n");
  std::cout << "Enter the number a: ";
  // printf("asdsasdsaasd3ewrfwe\n");
  std::cin >> input_a;
  // printf("sadsdasadsdas\n");

  // printf("d");

  // Declarando numeros
  if (mpz_set_str(N, input_N.c_str(), 10) != 0 || mpz_set_str(a, input_a.c_str(), 10) != 0) {
    std::cerr << "Invalid number format." << std::endl;
    mpz_clears(N, a, candidate, prime, generator, NULL);
    return 1;
  }

  // printf("e");

  // Comecando com um numero maior que N
  mpz_add_ui(candidate, N, 1);

  // printf("f");

  // Não contar numeros pares
  if (mpz_even_p(candidate)) {
    mpz_add_ui(candidate, candidate, 1);
  }

  // printf("g");


  while (true) {
    if (isPrimeMRtest(candidate, MR_count)) {
      mpz_set(prime, candidate);
      break;
    }
    mpz_add_ui(candidate, candidate, 2);
  }

  gmp_printf("Next prime: %Zd\n", prime);
  printf("Number of times the Miller-Rabin test was used: %d\n", MR_count);

  mpz_clears(N, a, candidate, prime, NULL);


  return 0;
}