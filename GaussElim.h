//
// Created by Yongjia Xu on 4/28/20.
//

#ifndef CS3353_PROJECT4_GAUSSELIM_H
#define CS3353_PROJECT4_GAUSSELIM_H

#include "Solver.h"
#include <cmath>
class GaussElim : public Solver{
private:
    double ** mat;
    int size;
public:
    GaussElim();
    bool forwardElim();
    void backSub();
    void solve();
    void print();
    void setMat(double**&, int);
};


#endif //CS3353_PROJECT4_GAUSSELIM_H
