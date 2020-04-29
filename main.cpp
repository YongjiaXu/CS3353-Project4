//
// Created by Yongjia Xu on 4/28/20.
//

#include <iostream>
#include "solverInterface.h"
using namespace std;

int main(){

    solverInterface s;
    s.setSolver(1);
    s.readIn("input.txt");
    s.solve();

    return 0;
}