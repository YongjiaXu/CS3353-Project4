//
// Created by Yongjia Xu on 4/28/20.
//

#include <sstream>
#include "solverInterface.h"
solverInterface::solverInterface(){
    this->mysolver = nullptr;
    this->size = 0;
}

bool solverInterface::readIn(string fileName) {
    ifstream file;
    file.open(fileName);
    if(!file.is_open()){
        cout << "Cannot open " << fileName << endl;
        return false;
    }
    string buffer;
    getline(file, buffer);
    size = stoi(buffer);
    mat = new double*[size];
    for(int i = 0; i < size; i++){
        mat[i] = new double [size+1];
    }
    int count = 0;
    while(count < size){
        if(file.eof()){
            break;
        }
        getline(file, buffer);
        if(buffer[0] == ' ' || buffer.empty()){
            break;
        }
        istringstream ss(buffer);
        vector<string> number;
        do {
            string n;
            ss >> n;
            number.push_back(n);
        } while (ss);
        for(int i = 0; i < size+1; i++){
            mat[count][i] = stod(number[i]);
        }
        count++;
    }
    for(int i = 0; i < size; i++){
        for(int z = 0; z < size+ 1; z++){
            cout << mat[i][z] << " ";
        }
        cout << endl;
    }
    file.close();
    return true;
}

void solverInterface::setSolver(int set) {
    if(set == 1){
        mysolver = new GaussElim;
    }
    else if(set == 2){
        mysolver == new GaussSeidel;
    }
    else if(set == 3){
        mysolver == new Jacobi;
    }
}

void solverInterface::clear() {
    if(mysolver != nullptr){
        delete mysolver;
    }
    for(int i = 0; i < size; i++){
        delete[] mat[i];
    }
    delete mat;
}