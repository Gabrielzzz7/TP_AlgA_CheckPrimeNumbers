#include <iostream>
#include <algorithm>
#include <gmpxx.h>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <random>
#include <ctime>

// ============= Implementação Miller-Rabin e achando menor primo maior que N ==============

bool PartialFactorization = false;

// Teste Miller-Rabin
// a^exp mod p
void expMod(mpz_class& result, const mpz_class& a, const mpz_class& exp, const mpz_class& primo) {
    mpz_class base = a % primo;
    mpz_class exp_copy = exp;
    result = 1;

    while (exp_copy > 0) {
        if (exp_copy % 2 == 1) { // if exp_copy is odd
            result = (result * base) % primo;
        }

        base = (base * base) % primo;
        exp_copy /= 2;
    }
}

mpz_class mod_exp(mpz_class base, mpz_class exp, mpz_class mod) {
    mpz_class result;
    expMod(result, base, exp, mod);
    return result;
}

// Teste Miller-Rabin
bool millerRabinTest(const mpz_class& n, const mpz_class& a) {
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
    expMod(x, a, d, n);

    // Se x == 1 ou x == n - 1, n é provavelmente primo
    if (x == 1 || x == n_minus_1) {
        return true;
    }

    // Realiza r - 1 iterações
    for (mpz_class i = 0; i < r - 1; i++) {
        // Calcula x = x^2 % n
        expMod(x, x, 2, n);

        // Se x == n - 1, n é provavelmente primo
        if (x == n_minus_1) {
            return true;
        }

        // Se x == 1, n é composto
        if (x == 1) {
            return false;
        }
    }

    return false; // n é composto
}

// Teste de primalidade de Miller-Rabin
bool isPrimeMRtest(const mpz_class& n, int& MR_count) {
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

std::pair<mpz_class, mpz_class> findEstimativeGenerator(mpz_class prime, std::vector<std::vector<mpz_class>> factors) {
    mpz_class n, result, exp, q, b, ord;
    n = prime - 1;
    b = 1;
    ord = 1;
    for (const auto& factor_group : factors) {
        for (int i = 2; i < n; i++) {
            q = factor_group[0];
            exp = n / q;
            expMod(result, i, exp, prime);
            if (result != 1) {
                b *= result;
                ord *= q;
                break;
            }
        }
    }

    b %= prime;
    ord %= prime;
    return { b, ord };
}

mpz_class findPrimitiveRoot(mpz_class primo, std::vector<std::vector<mpz_class>> factors) {
    mpz_class n, result, exp, q;

    n = primo - 1;

    for (mpz_class gerador = 1; gerador < primo; gerador++) {
        bool is_primitive = true;

        for (const auto& factor_group : factors) {

            q = factor_group[0];
            exp = n / q;
            expMod(result, gerador, exp, primo);

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

// ========== FORÇA BRUTA =============
std::pair<mpz_class, double> bruteForce(mpz_class p, mpz_class g, mpz_class a) {
    clock_t begin = clock();
    const int limitTime = 120;

    mpz_class i;
    mpz_class pot, modulo;

    pot = 1; // g ^ 0

    for (i = 0; i < p; i++) {
        if ((clock() - begin) / CLOCKS_PER_SEC > limitTime) {
            return std::make_pair(-1, (double)(clock() - begin) / CLOCKS_PER_SEC);
        }
        mpz_mod(modulo.get_mpz_t(), pot.get_mpz_t(), p.get_mpz_t());

        if (modulo == a) {
            clock_t end_time = clock();
            double elapsed_time = double(end_time - begin) / CLOCKS_PER_SEC; // Tempo decorrido em segundos

            if (elapsed_time > limitTime) {
                return { -1, elapsed_time };
            }

            return { i, elapsed_time };
        }

        pot *= g;

        if (pot > p) {
            mpz_mod(pot.get_mpz_t(), pot.get_mpz_t(), p.get_mpz_t());
        }
    }

    clock_t end_time = clock(); // Tempo de término
    double elapsed_time = double(end_time - begin) / CLOCKS_PER_SEC;

    if (elapsed_time > limitTime) {
        return { -1, elapsed_time };
    }

    return { -1, elapsed_time };
}

// ========== BSGS =============
std::pair<mpz_class, double> BSGS(const mpz_class& p, const mpz_class& g, const mpz_class& a) {
    // Obtendo o tempo de início do processamento 
    clock_t begin = clock();
    const int limitTime = 120;
    mpz_class r; // Teto da raiz do número primo p
    mpz_class c; // Resto da divisão de g ^ r por p
    mpz_class fat1, fat2; // factors que iremos comparar nas iterações
    mpz_class u, l; // Expoentes dos factors
    mpz_class ans;

    std::vector<mpz_class> factors;
    std::vector<mpz_class>::iterator it;

    // Calculando a raiz de p
    mpz_sqrt(r.get_mpz_t(), p.get_mpz_t());
    r = r + 1; // Definindo o teto

    // Calculando g^r % p
    expMod(c, g, r, p);

    // Pré-cálculo de todos os valores de g^u * a % p para 0 <= u < r
    for (u = 0; u < r; u++) {
        // Calculando a * g^u % p
        expMod(fat2, g, u, p);
        fat2 *= a;
        mpz_mod(fat2.get_mpz_t(), fat2.get_mpz_t(), p.get_mpz_t());

        // Verificando se o tempo de processamento excedeu o limite
        if ((clock() - begin) / CLOCKS_PER_SEC > limitTime) {
            return std::make_pair(-1, (double)(clock() - begin) / CLOCKS_PER_SEC);
        }

        factors.push_back(fat2);
    }

    expMod(fat1, c, r, p);

    // Buscando por colisão entre c^u % p e os valores de a * g^u % p pré-calculados
    for (u = 0; u < r; u++) {
        if ((clock() - begin) / CLOCKS_PER_SEC > limitTime) {
            return std::make_pair(-1, (double)(clock() - begin) / CLOCKS_PER_SEC);
        }

        expMod(fat1, c, u, p);
        it = std::find(factors.begin(), factors.end(), fat1);

        if (it != factors.end()) {
            l = it - factors.begin();
            ans = u * r - l;
            if (ans < p)
                return std::make_pair(ans, (double)(clock() - begin) / CLOCKS_PER_SEC);
        }
    }

    return std::make_pair(-1, (double)(clock() - begin) / CLOCKS_PER_SEC);
}

// ========== TEOREMA RESTO CHINES =========

// Função para calcular o inverso modular
mpz_class modInverse(const mpz_class& a, const mpz_class& m) {
    mpz_class inv;
    mpz_invert(inv.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t());
    return inv;
}

// Função que verifica se dois números são coprimos
bool areCoprime(const mpz_class& a, const mpz_class& b) {
    mpz_class gcd;
    mpz_gcd(gcd.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return gcd == 1;
}

// Função que implementa o Teorema do Resto Chinês

mpz_class chineseRest(const std::vector<std::vector<mpz_class>>& n, const std::vector<mpz_class>& a) {
    // Calcula o produto de todos os módulos
    mpz_class prod = 1;
    for (const auto& vec : n) {
        mpz_class prime = vec[0];
        mpz_pow_ui(prime.get_mpz_t(), prime.get_mpz_t(), vec.size());
        prod *= prime;
    }

    // Aplica o Teorema do Resto Chinês
    mpz_class result = 0;
    for (size_t i = 0; i < n.size(); ++i) {
        mpz_class ni = n[i][0];
        mpz_class ai = a[i];
        mpz_class p = prod / ni;
        mpz_class inv = modInverse(p, ni);
        result += ai * inv * p;
    }

    result %= prod;
    return result;
}

// ========== POHLIG-HELLMAN ===========

std::pair<mpz_class, double> pollingHellman(const mpz_class& a, const mpz_class& p, const mpz_class& g, std::vector<std::vector<mpz_class>> factors) {
    clock_t begin = clock();
    const int limitTime = 120;
    std::vector<mpz_class> xValues;
    mpz_class n = p - 1;

    for (const auto& qFactors : factors) {
        if ((clock() - begin) / CLOCKS_PER_SEC > limitTime) {
            return std::make_pair(-1, (double)(clock() - begin) / CLOCKS_PER_SEC);
        }
        mpz_class pi = qFactors[0];
        mpz_class ei = qFactors.size();

        mpz_class gi;
        expMod(gi, g, n / pi, n);

        // Passo 2: Calcular h_i
        mpz_class hi;
        expMod(hi, p, n / pi, n);

        // Passo 3: Resolver a congruência g_i^x_i ≡ h_i (mod p_i^e_i)
        mpz_class xi;
        expMod(xi, hi, (pi - 1) * (ei - 1), n);

        xValues.push_back(xi);
    }
    mpz_class result = chineseRest(factors, xValues);

    return std::make_pair(result, (double)(clock() - begin) / CLOCKS_PER_SEC);
}

// ========== Main =============
int main() {
    mpz_class N, a, candidate, prime, generator, discrete_log;
    int MR_count = 0;
    double time;
    std::string input_N, input_a;
    std::chrono::seconds time_limit = std::chrono::minutes(1);

    std::cout << "Digite um numero grande N: ";
    std::getline(std::cin, input_N);

    std::cout << "Insira o valor para a: ";
    std::getline(std::cin, input_a);

    // Converter inputs para mpz_class
    N = input_N;
    a = input_a;

    // Começando com um número maior que N
    candidate = N + 1;

    // Não contar números pares
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
        generator = findPrimitiveRoot(prime, factors);
        std::cout << "Gerador g de " << prime << " é: ";
        mpz_out_str(stdout, 10, generator.get_mpz_t());
        std::cout << std::endl;
    }

    //Logaritmo discreto
    std::pair<mpz_class, double> result_pair = pollingHellman(a, prime, generator, factors);
    discrete_log = result_pair.first;
    time = result_pair.second;

    if (discrete_log == -1) {
        std::cout << "O tempo de execução execedeu 2 minutos" << std::endl;
    }
    else {
        std::cout << "O logaritmo discreto de " << a << " módulo " << prime << " na base " << generator << " é: ";
        mpz_out_str(stdout, 10, discrete_log.get_mpz_t());
        std::cout << std::endl;
        std::cout << "O tempo de execução foi: " << time << " segundos" << std::endl;
    }

    return 0;
}