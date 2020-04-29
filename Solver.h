//
// Created by Yongjia Xu on 4/28/20.
//

#ifndef CS3353_PROJECT4_SOLVER_H
#define CS3353_PROJECT4_SOLVER_H

#include <iostream>
using namespace std;
class Solver {
public:
    Solver(){}
    virtual void solve() = 0;
    virtual void setMat(double**&, int) = 0;
};


#endif //CS3353_PROJECT4_SOLVER_H
