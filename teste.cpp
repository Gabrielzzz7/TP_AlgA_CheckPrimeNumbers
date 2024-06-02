#include <iostream>
#include <algorithm>
#include <gmpxx.h>
#include <vector>
#include <random>
#include <ctime>

// ============= Implementação Miller-Rabin e achando menor primo maior que N ==============

bool PartialFactorization = false;

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

std::vector<std::vector<mpz_class>> factorize(mpz_class n) {
    std::vector<std::vector<mpz_class>> factorization;

    mpz_class d;
    for (d = 2; d * d <= n && d < 1e7; d++) {
        while (n % d == 0) {
            bool found = false;
            for (auto& factor_group : factorization) {
                if (factor_group[0] == d) {
                    factor_group.push_back(d);
                    found = true;
                    break;
                }
            }
            if (!found) {
                factorization.push_back({ d });
            }
            n /= d;
        }
    }
    if (d >= 1e7) {
        PartialFactorization = true;
    }
    if (n > 1) {
        bool found = false;
        for (auto& factor_group : factorization) {
            if (factor_group[0] == n) {
                factor_group.push_back(n);
                found = true;
                break;
            }
        }
        if (!found) {
            factorization.push_back({ n });
        }
    }
    return factorization;
}

mpz_class randInt(mpz_class min, mpz_class max) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<mpz_class> dis(min.get_ui(), max.get_ui());
    return mpz_class(dis(gen));
}

std::pair<mpz_class, mpz_class> findEstimativeGenerator(mpz_class prime, std::vector<std::vector<mpz_class>> factors) {
    mpz_class n, b, ord;
    n = prime - 1;
    b = 1;
    ord = 1;


}

mpz_class find_primitive_root(mpz_class primo, std::vector<std::vector<mpz_class>> factors) {
    std::vector<std::vector<mpz_class>> factors;
    mpz_class n, result, exp, q;

    n = primo - 1;

    for (mpz_class gerador = 1; gerador < primo; gerador++) {

        bool is_primitive = true;

        for (const auto& factor_group : factors) {

            q = factor_group[0];
            exp = n / q;
            mpz_powm(result.get_mpz_t(), gerador.get_mpz_t(), exp.get_mpz_t(), primo.get_mpz_t());

            if (result == 1) {
                is_primitive = false;
                break;
            }
        }

        if (is_primitive) {
            // Se chegou aqui, gerador é primitivo
            return gerador;
        }
    }

    return -1;
}

// ========== BSGS =============

mpz_class BSGS(mpz_class p, mpz_class g, mpz_class a) {

    mpz_class r; // Teto da raiz do número primo p
    mpz_class c; // Resto da divisão de g ^ r por p
    mpz_class fat1, fat2; // Fatores que iremos comparar nas iterações
    mpz_class u, l; // Expoentes dos fatores
    mpz_class ans;

    std::vector<mpz_class> fatores;
    std::vector<mpz_class>::iterator it;

    mpz_sqrt(r.get_mpz_t(), p.get_mpz_t()); // Calculando raiz de p
    r = r + 1; // Definindo o teto
    mpz_powm(c.get_mpz_t(), g.get_mpz_t(), r.get_mpz_t(), p.get_mpz_t());

    for (u = 0; u < r; u++) {

        // a * g ^ u mod p
        mpz_pow_ui(fat2.get_mpz_t(), g.get_mpz_t(), (unsigned long)u.get_mpz_t());
        fat2 *= a;
        mpz_mod(fat2.get_mpz_t(), fat2.get_mpz_t(), p.get_mpz_t());

        fatores.push_back(fat2);

        // c ^ l mod p
        mpz_powm(fat1.get_mpz_t(), c.get_mpz_t(), u.get_mpz_t(), p.get_mpz_t());
        it = std::find(fatores.begin(), fatores.end(), fat2);

        if (it != fatores.end()) {

            l = (mpz_class)(it - fatores.begin());
            ans = l * r - u;
            return ans;
        }
    }

    return -1;
}


// ========== Main =============

int main() {
    mpz_class N, a, candidate, prime, generator;
    int MR_count = 0;

    // Leitura
    std::string input_N, input_a;

    std::cout << "Digite um numero grande N: ";
    std::getline(std::cin, input_N);

    std::cout << "Insira o valor para a: ";
    std::getline(std::cin, input_a);

    // Convert inputs to mpz_class
    N = input_N;
    a = input_a;

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

    std::vector<std::vector<mpz_class>> factors;

    factors = factorize(prime - 1);

    if (PartialFactorization) {
        mpz_class minOrder;
        std::pair<mpz_class, mpz_class> estimative = findEstimativeGenerator(prime, factors);
        generator = estimative.first;
        minOrder = estimative.second;
        std::cout << "Estimativa do gerador: ";
        mpz_out_str(stdout, 10, generator.get_mpz_t());
        std::cout << "\nOrdem mínima estimada: ";
        mpz_out_str(stdout, 10, minOrder.get_mpz_t());
        std::cout << std::endl;

    }
    else {
        // Raiz primitiva
        generator = find_primitive_root(prime, factors);
    }

    std::cout << "Gerador g de " << prime << " e: ";
    mpz_out_str(stdout, 10, generator.get_mpz_t());
    std::cout << std::endl;

    return 0;
}