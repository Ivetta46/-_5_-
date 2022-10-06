#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <iostream>
#include <utility>
#include <vector>

class Borders{
public:
    double min;
    double max;

    Borders(){};
    Borders(double min, double max);
};


class Polynomial{
private:
    double* coef;
    double** d_coef;
    int N;
    void makeDerivative();

public:
    std::vector<double> roots;
    std::vector<Borders> segments;
    double epsilon = 0.001;
    Polynomial();
    Polynomial(int n, double* coef);

    Borders getRing();
    int getDegree();

    double derivative(int order, double arg);
    void printDeriv(int order);

    void dihotomia(double left, double right);
    int isSingleRoot(double left, double right);
    void findRootPrecise(double left, double right);

    void lokalizeRoots();
    void findRoots();
    void printRoots();

    ~Polynomial();
};



#endif