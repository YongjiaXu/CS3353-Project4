//
// Created by Yongjia Xu on 4/28/20.
//

#include <iostream>
#include "solverInterface.h"
using namespace std;

int main(){

    solverInterface s;
    s.setSolver(3);
    s.readIn("input.txt");

    return 0;
}