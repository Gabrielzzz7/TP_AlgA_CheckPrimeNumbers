#include <iostream>
#include <algorithm>
#include <gmpxx.h>
#include <vector>

// ============= Implementação Miller-Rabin e achando menor primo maior que N ==============

// a^exp mod p
void exp_mod(mpz_class& result, const mpz_class& a, const mpz_class& exp, const mpz_class& primo) {
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
    exp_mod(x, a, d, n);

    // Se x == 1 ou x == n - 1, n é provavelmente primo
    if (x == 1 || x == n_minus_1) {
        return true;
    }

    // Realiza r - 1 iterações
    for (mpz_class i = 0; i < r - 1; i++) {
        // Calcula x = x^2 % n
        exp_mod(x, x, 2, n);

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

mpz_class find_primitive_root(const mpz_class& primo) {
    std::vector<mpz_class> factors;
    mpz_class n, result, exp, q;

    n = primo - 1;

    // Fatorando primo-1
    factors = factorize(n);

    for (mpz_class gerador = 1; gerador < primo; gerador++) {
        bool is_primitive = true;

        for (size_t i = 0; i < factors.size(); i++) {
            q = factors[i];
            exp = n / q;
            exp_mod(result, gerador, exp, primo);

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

mpz_class brute_force(mpz_class p, mpz_class g, mpz_class a){

    mpz_class i;
    mpz_class pot, modulo;

    pot = 1; // g ^ 0

    for (i = 0; i < p; i++){
        
        mpz_mod(modulo.get_mpz_t(), pot.get_mpz_t(), p.get_mpz_t());
        
        if(modulo == a){
            return i;
        }
        
        pot *= g;

        if(pot > p){

            mpz_mod(pot.get_mpz_t(), pot.get_mpz_t(), p.get_mpz_t()); 
        }        
    }
    
    return -1;
}

// ========== BSGS =============

mpz_class BSGS(const mpz_class& p, const mpz_class& g, const mpz_class& a) {
    
    mpz_class r; // Teto da raiz do número primo p
    mpz_class c; // Resto da divisão de g ^ r por p
    mpz_class fat1, fat2; // Fatores que iremos comparar nas iterações
    mpz_class u, l; // Expoentes dos fatores
    mpz_class ans;

    std::vector<mpz_class> fatores;
    std::vector<mpz_class>::iterator it;

    mpz_sqrt(r.get_mpz_t(), p.get_mpz_t()); // Calculando raiz de p
    r = r + 1; // Definindo o teto
    exp_mod(c, g, r, p);

    for (u = 0; u < r; u++) {
        // a * g ^ u mod p
        mpz_pow_ui(fat2.get_mpz_t(), g.get_mpz_t(), u.get_ui()); // g ^ u
        fat2 *= a; // Multiplicando por a
        mpz_mod(fat2.get_mpz_t(), fat2.get_mpz_t(), p.get_mpz_t()); // Calculando resto da divisão por p

        fatores.push_back(fat2);

        // c ^ l mod p
        exp_mod(fat1, c, u, p);
        it = std::find(fatores.begin(), fatores.end(), fat2);

        if (it != fatores.end()) {
            l = it - fatores.begin();
            ans = l * r - u;
            return ans;
        }
    }

    return -1;
}

// ========== Main =============
int main() {

    mpz_class N, a, candidate, prime, generator, discrete_log;
    int MR_count = 0;
    std::string input_N, input_a;

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

    // Raiz primitiva
    generator = find_primitive_root(prime);

    std::cout << "Gerador g de " << prime << " é: ";
    mpz_out_str(stdout, 10, generator.get_mpz_t());
    std::cout << std::endl;

    //Logaritmo discreto
    discrete_log = brute_force(prime, generator, a);

    std::cout << "O logaritmo discreto de " << a << " módulo " << prime << " na base " << generator << " é: ";
    mpz_out_str(stdout, 10, discrete_log.get_mpz_t());
    std::cout << std::endl;

    return 0;
}