#include "polynomial.h"
#include <iomanip>
#include <math.h>


void Polynomial::dihotomia(double left, double right){

    if (right - left > this->epsilon && isSingleRoot(left, right) > 1 ){
        dihotomia(left, (left + right) / 2.);
        dihotomia((left + right) / 2., right);
    }

    if (isSingleRoot(left, right) == 1){
        this->segments.push_back({left, right});
    }

}


int Polynomial::isSingleRoot(double left, double right){

    int count1 = 0;
    int count2 = 0;

    for (int i = 1; i <= this->N; i++){

        if (signbit(this->derivative(i - 1, left) * this->derivative(i, left))){
            count1++;
        }

        if (signbit(this->derivative(i - 1, right) * this->derivative(i, right))){
            count2++;
        }

    }

    return count1 - count2;

}

void Polynomial::lokalizeRoots(){
    this->dihotomia(getRing().min, getRing().max);
}


void Polynomial::findRootPrecise(double left, double right){
    if (right - left >= this->epsilon && 100*this->derivative(0, left) * this->derivative(0, right) < 0){
        this->findRootPrecise(left, (left + right) / 2.);
        this->findRootPrecise((left + right) / 2., right);
    }
    if(right - left < this->epsilon && 100*this->derivative(0, left) * this->derivative(0, right) < 0){
        this->roots.push_back((left + right) / 2.);
    }

}


void Polynomial::findRoots(){
    for (auto i: this->segments){
        this->findRootPrecise(i.min, i.max);
    }
}


void Polynomial::printRoots(){

    double tmp = this->epsilon;
    int count = 0;

    while (tmp < 1){
        tmp *= 10;
        count++;
    }

    std::cout << std::setprecision(count) << std::fixed;
    for (auto i: this->roots){
        std::cout << i << " ";
    }
    std::cout << std::endl;

}

