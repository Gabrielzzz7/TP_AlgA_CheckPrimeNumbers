#include <iostream>
#include <gmpxx.h>
#include <vector>
// ============= Implementação Miller-Rabin e achando menor primo maior que N ==============

// Teste Miller-Rabin
bool millerRabinTest(mpz_class n, mpz_class a) {
    mpz_class x, n_minus_1, d, r;

    n_minus_1 = n - 1;

    // Inicializa d como n - 1
    d = n_minus_1;

    // Calcula o valor de d e r tal que d * 2^r = n - 1
    r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r += 1;
    }

    // Calcula x = a^d % n
    mpz_powm(x.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

    // Se x == 1 ou x == n - 1, n é provavelmente primo
    if (x == 1 || x == n_minus_1) {
        return true;
    }

    // Realiza r - 1 iterações
    mpz_class i;
    i = 0;
    while (i < r) {
        // Calcula x = x^2 % n
        mpz_powm_ui(x.get_mpz_t(), x.get_mpz_t(), 2, n.get_mpz_t());

        // Se x == n - 1, n é provavelmente primo
        if (x == n_minus_1) {
            return true;
        }

        // Se x == 1, n é composto
        if (x == 1) {
            return false;
        }

        i += 1;
    }

    return false; // n é composto
}

// Teste de primalidade de Miller-Rabin
bool isPrimeMRtest(mpz_class n, int& MR_count) {
    mpz_class primes[40];

    // Inicializa os números primos
    const char* primes_str[40] = {
        "2", "3", "5", "7", "11", "13", "17", "19", "23", "29",
        "31", "37", "41", "43", "47", "53", "59", "61", "67", "71",
        "73", "79", "83", "89", "97", "101", "103", "107", "109", "113",
        "127", "131", "137", "139", "149", "151", "157", "163", "167", "173"
    };

    for (int i = 0; i < 40; ++i) {
        primes[i] = mpz_class(primes_str[i]);
    }

    bool isPrime = true;

    for (int i = 0; i < 40; i++) {
        MR_count++;
        if (!millerRabinTest(n, primes[i])) { // Usa os 40 primeiros primos
            isPrime = false;
            break;
        }
    }

    return isPrime;
}

// =========== Implementação - achando raíz primitiva de N ============

std::vector<mpz_class> factorize(mpz_class n) {
    
    std::vector<mpz_class> factorization;

    for (mpz_class d = 2; d * d <= n; d++) {
        while (n % d == 0) {
            factorization.push_back(d);
            n /= d;
        }
    }
    if (n > 1)
        factorization.push_back(n);
    return factorization;
}

void find_primitive_root(mpz_class& gerador, mpz_class& primo) {

    std::vector<mpz_class> factors;
    mpz_class n, result, exp;

    n = primo - 1;

    // Fatorando primo-1
    factors = factorize(n);
    for (gerador = 1; gerador < primo; ++gerador) {
        bool is_primitive = true;
        for (size_t i = 0; i < factors.size(); ++i) {
            mpz_class q = factors[i];
            exp = n / q;
            mpz_powm(result.get_mpz_t(), gerador.get_mpz_t(), exp.get_mpz_t(), primo.get_mpz_t());
            if (result == 1) {
                is_primitive = false;
                break;
            }
        }
        if (is_primitive) {
            // Se chegou aqui, gerador é primitivo
            break;
        }
    }
}



// ========== Main =============

int main() {
    mpz_class N, a, candidate, prime, generator;
    int MR_count = 0;

    // Leitura
    std::string input_N;
    std::cout << "Digite um numero grande N: ";
    std::getline(std::cin, input_N);
    
    // Convert input to mpz_class
    N = input_N;

    // Comecando com um numero maior que N
    candidate = N + 1;

    // Não contar numeros pares
    if (candidate % 2 == 0) {
        candidate += 1;
    }

    while (true) {
        if (isPrimeMRtest(candidate, MR_count)) {
            prime = candidate;
            break;
        }
        candidate += 2;
    }

    gmp_printf("Proximo primo: %Zd\n", prime.get_mpz_t());
    printf("Numero de testes de Miller-Rabin: %d\n", MR_count);

    // Raiz primitiva
    find_primitive_root(generator, prime);

    std::cout << "Gerador g de " << prime << " e: ";
    mpz_out_str(stdout, 10, generator.get_mpz_t());
    std::cout << std::endl;

    return 0;
}