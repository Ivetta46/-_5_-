#include "polynomial.h"
#include <cmath>

Polynomial::Polynomial(){};

Polynomial::Polynomial(int n, double* coef){
    this->N = n;
    this->coef = coef;
    this->makeDerivative();
}

int Polynomial::getDegree(){
    return this->N;
}

Borders Polynomial:: getRing(){

    double A = -1;
    double B = -1;

    for (int i = 0; i <= this->N; i++){

        if (A < sqrt(pow(coef[i], 2)) && i > 0){
            A = sqrt(pow(coef[i], 2));
        }

        if (B < sqrt(pow(coef[i], 2)) && i < this->N){
            B = sqrt(pow(coef[i], 2));
        }
    }

    double max, min;

    min = sqrt(pow(coef[this->N], 2)) / (sqrt(pow(coef[this->N],2)) + B);
    max = 1.0 + A / sqrt(pow(coef[0], 2));

    std::cout << min << " " << max << std::endl;
    return {min, max};

}

Borders::Borders(double min, double max){
    this->min = min;
    this->max = max;
}


double Polynomial::derivative(int order, double arg){

    double deriv = 0;
    double tmp = 1;

    for (int i = 0; i <= this->N; i++){

        tmp = this->d_coef[order][i] * pow(arg, this->N - i);

        /*	for (int j = 0; j < this->N - i; j++){
                tmp *= arg;
            }
    */
        deriv += tmp;
    }


    return deriv;

}


void Polynomial::makeDerivative(){
    this->d_coef = new double*[this->N + 1];

    for (int i = 0; i <= this->N; i++){
        this->d_coef[i] = new double[this->N + 1];
    }

    for (int i = 0; i <= this->N; i++){
        this->d_coef[0][i] = this->coef[i];
    }

    for (int i = 1; i <= this->N; i++){
        for (int j = 0; j <= this->N; j++){
            this->d_coef[i][j] = this->d_coef[i - 1][j] * (this->N - j);
        }

        for (int j = this->N - 1; j >= 0; j--){
            this->d_coef[i][j + 1] = this->d_coef[i][j];
        }

        this->d_coef[i][0] = 0;
    }
}


Polynomial::~Polynomial(){
    //delete[] this->coef;

    for (int i = 0; i <= this->N; i++){
        delete[] this->d_coef[i];
    }
}


void Polynomial::printDeriv(int order){

    for (int i = 0; i <= this->N; i++){
        std::cout << this->d_coef[order][i] << " ";
    }

    std::cout << std::endl;
}