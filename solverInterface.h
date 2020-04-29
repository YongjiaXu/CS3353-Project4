//
// Created by Yongjia Xu on 4/28/20.
//

#ifndef CS3353_PROJECT4_SOLVERINTERFACE_H
#define CS3353_PROJECT4_SOLVERINTERFACE_H

#include "Solver.h"
#include "GaussElim.h"
#include "GaussSeidel.h"
#include "Jacobi.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
class solverInterface {
private:
    Solver *mysolver;
    double **mat;
    int size;
public:
    solverInterface();
    bool readIn(string);
    void setSolver(int);
    void solve();
    void clear();

};


#endif //CS3353_PROJECT4_SOLVERINTERFACE_H
