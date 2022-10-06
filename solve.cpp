#include <iostream>
#include "polynomial.h"
#include <cmath>
#include <iomanip>
#include <clocale>
#include <windows.h>


int main(){
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    double gamma0, ro0, P0, U0;
    double gamma3, C3, P3, U3;

    std::cout << "Enter precision epsilon." << std::endl;

    double epsilon;
    std::cin >> epsilon;

    std::cout << std::endl;

    std::cout << "Enter gamma, rho, P, U  (0)" << std::endl;

    std::cin >> gamma0;
    std::cin >> ro0;
    std::cin >> P0;
    std::cin >> U0;

    std::cout << std::endl;

    std::cout << "Enter gamma, C, P, U (3)" << std::endl;

    std::cin >> gamma3;
    std::cin >> C3;
    std::cin >> P3;
    std::cin >> U3;

    std::cout << std::endl;



    double e0 = P0/(ro0 * (gamma0 - 1));
    double ro3 = (P3 * gamma3)/(C3 * C3);
    double X = P3 / P0;
    double n = 2 * gamma3 / (gamma3 - 1);
//    double Z = pow((P1/P3), (1/n));
    double alpha0 = (gamma0 + 1) / (gamma0 - 1);
    double mu = (U3 - U0) * sqrt(((gamma0 - 1) * ro0)/2 / P0);

    double nu = (2/(gamma3 - 1) * sqrt((gamma3 * (gamma0 - 1)/2) * (P3/P0) * (ro0/ro3)));



    double *coef = new double[2 * (int)n];
    memset(coef, 0, (2 * n) * sizeof(double));
    coef[0] = pow(X, 2);
    coef[(int)n - 2] = -alpha0 * nu * nu * X;
    coef[(int)n - 1] = 2 * alpha0 * nu * (nu + mu) * X;
    coef[(int)n] = -(2 + pow((mu + nu), 2) * alpha0) * X;
    coef[2 * (int)n - 2] = -nu * nu;
    coef[2 * (int)n - 1] = 2 * nu * (mu + nu);
    coef[2 * (int)n] = -pow((mu + nu), 2) + 1;

    std::cout << "2n coef: " << coef[2 * (int)n] << std::endl;

    Polynomial pol = Polynomial(2 * n, coef);
    pol.epsilon = epsilon;
    pol.lokalizeRoots();
    pol.findRoots();


    std::cout << std::setprecision(8) << std::fixed;
    std::cout << "Coefs:"<< std::endl;
    int count = 0;
    for (int i = 0; i < 2 * n; i++){
        if (coef[i] != 0){
            std::cout << count << " " << coef[i] << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Intervals:" << std::endl;
    for (auto i: pol.segments){
        std::cout << "(" << i.min << ",   " << i.max << ")";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;


    std::cout << "Roots:" << std::endl;
    pol.printRoots();
    std::cout << std::endl;

    std::cout << "Speeds:" << std::endl;
    for (auto root: pol.roots){
        std::cout << "Roots: " << root << std::endl;
        double P1 = P3 * pow(root, n);
        double s = P1 / P0;
        double C0 = pow(gamma0 * (P0 / ro0), 0.5);
        double U_plus_1 = U0 + ((C0 * ((P1/P3) * (P3/P0) - 1))/pow(0.5 * gamma0 * (gamma0 - 1) * (1 + ((gamma0 + 1)/(gamma0 - 1)) * (P1/P3) * (P3/P0)), 0.5));
        double U_minus_1 = U0 - ((C0 * ((P1/P3) * (P3/P0) - 1))/pow(0.5 * gamma0 * (gamma0 - 1) * (1 + ((gamma0 + 1)/(gamma0 - 1)) * (P1/P3) * (P3/P0)), 0.5));
        std::cout << "U1 = " << U_plus_1 << " ||  " << "U1 = " << U_minus_1 << std::endl;

        double rho_frac = ((gamma0 - 1) + (gamma0 + 1) * s) / ((gamma0 + 1) + (gamma0 - 1) * s); //ro1/ro0
        double D0_plus = U0 + sqrt(rho_frac * (P1 - P0) / (rho_frac * ro0 - ro0));
        double D0_minus = U0 - sqrt(rho_frac * (P1- P0) / (rho_frac * ro0 - ro0));
        std::cout << "D0 = " << D0_plus << " ||  " << "D0 = " << D0_minus << std::endl;
    }

    return 0;
}