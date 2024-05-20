#include <iostream>
#include <gmp.h>

// a principio, o codigo funciona para grande parte dos problemas pedidos. A ultima parte ficou problematica (3 questao)

// chcando se o numero e primo
bool is_prime(const mpz_t n, int iterations) {
    return mpz_probab_prime_p(n, iterations) > 0;
}

// achando maior divisor comum
void mdc(mpz_t result, const mpz_t a, const mpz_t b) {
    mpz_mdc(result, a, b);
}

// checando se g e raiz primitiva modulo p
bool is_primitive_root(const mpz_t g, const mpz_t p) {
    mpz_t order, exp, gcd_result, temp;
    mpz_inits(order, exp, gcd_result, temp, NULL);
    
    // p - 1
    mpz_sub_ui(order, p, 1);
    
    // Check for each factor of p-1
    mpz_tdiv_q_ui(exp, order, 2);
    mdc(mdc_result, g, p);
    
    if (mpz_cmp_ui(mdc_result, 1) != 0) {
        mpz_clears(order, exp, mdc_result, temp, NULL);
        return false;
    }
    
    mpz_powm(temp, g, exp, p);
    if (mpz_cmp_ui(temp, 1) == 0) {
        mpz_clears(order, exp, mdc_result, temp, NULL);
        return false;
    }
    
    mpz_clears(order, exp, mdc_result, temp, NULL);
    return true;
}

// Funcao raiz primitiva modulo p
void find_primitive_root(mpz_t result, const mpz_t p) {
    mpz_t g;
    mpz_init_set_ui(g, 2); // Primeiro candidato para g
    
    while (true) {
        if (is_primitive_root(g, p)) {
            mpz_set(result, g);
            break;
        }
        mpz_add_ui(g, g, 1);
    }
    
    mpz_clear(g);
}

int main() {
    mpz_t N, a, candidate, prime, generator;
    int iterations = 25; // valor de confianca
    int miller_rabin_count = 0;

    mpz_init(N);
    mpz_init(a);
    mpz_init(candidate);
    mpz_init(prime);
    mpz_init(generator);

    // Leitura
    std::string input_N, input_a;
    std::cout << "Enter a large number N: ";
    std::cin >> input_N;
    std::cout << "Enter the number a: ";
    std::cin >> input_a;

    // Declarando numeros
    if (mpz_set_str(N, input_N.c_str(), 10) != 0 || mpz_set_str(a, input_a.c_str(), 10) != 0) {
        std::cerr << "Invalid number format." << std::endl;
        mpz_clears(N, a, candidate, prime, generator, NULL);
        return 1;
    }

    // Comecando com um numero maior que N
    mpz_add_ui(candidate, N, 1);

    // Não contar numeros pares
    if (mpz_even_p(candidate)) {
        mpz_add_ui(candidate, candidate, 1);
    }

    // Achar menor primo maior que N
    while (!is_prime(candidate, iterations)) {
        mpz_add_ui(candidate, candidate, 2);
        miller_rabin_count++;
    }
    
    mpz_set(prime, candidate);

    // Raiz primitiva 
    find_primitive_root(generator, prime);

    // Output
    std::cout << "Menor primo maior que " << input_N << " e: ";
    mpz_out_str(stdout, 10, prime);
    std::cout << std::endl;

    std::cout << "Raiz primitiva modulo " << input_N << " e: ";
    mpz_out_str(stdout, 10, generator);
    std::cout << std::endl;

    std::cout << "Miller-Rabin foi usado" << miller_rabin_count << " vezes." << std::endl;

    // clear
    mpz_clears(N, a, candidate, prime, generator, NULL);

    return 0;
}