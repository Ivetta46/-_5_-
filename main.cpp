#include <iostream>
#include "float.h"
#include <cmath>
#include <utility>
#include "polynomial.h"


int main(){
    int n;

    std::cin >> n;

    double* coef = new double[n+1];


    for (int i = 0; i <= n; i++){
        std::cin >> coef[i];
    }

    Polynomial pol = Polynomial(n, coef);

    Borders borders = pol.getRing();

    std::cout << borders.min << " " << borders.max << std::endl;

    pol.lokalizeRoots();
    pol.findRoots();

    pol.printRoots();


    return 0;
}
