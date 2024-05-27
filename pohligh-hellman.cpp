#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<long long> fatorar(long long n) {

    vector<long long> factorization;

    for (long long d = 2; d * d <= n; d++) {
        while (n % d == 0) {
            factorization.push_back(d);
            n /= d;
        }
    }
    if (n > 1)
        factorization.push_back(n);
    return factorization;
}

// a --> logaritmando; g --> base do logaritmo
long long BSGS(long long p, long long g, long long a){

    long long r; // Teto da raiz do número primo p
    long long c; // Resto da divisão de g ^ r por p
    long long fat1, fat2; // Fatores que iremos comparar nas iterações

    r = (long long) ceil(sqrtl(p));
    c = (long long) powl(g, r);

    for(int i = 0; i < r; i++){

        fat1 = (long long) powl(c, i);
        fat2 = a * ( (long long) powl(g, i) );

        if(fat1 == fat2){
            break;
        }
    }

}