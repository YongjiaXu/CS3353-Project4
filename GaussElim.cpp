//
// Created by Yongjia Xu on 4/28/20.
//

#include "GaussElim.h"
GaussElim::GaussElim(){
    this->mat = nullptr;
    this->size = 0;
}

void GaussElim::setMat(double**& mat, int size) {
    this->mat = mat;
    this->size = size;
}

void GaussElim::solve() {
    //gauss elimination: decomposition to upper triangle matrix
    //also check if the matrix is singular - non-invertible thus have no unique solution
    bool singular = forwardElim();
    // if the matrix is singular
    if (singular){
        cout << "This is a singular matrix.\n" << endl;

        if (mat[singular][size]){
            cout << "The matrix is a singular matrix." << endl;
        }
            //right hand size is all zero
        else{
            cout << "May have infinitely many solutions." << endl;
        }
        return;
    }

    //backward substitution
    backSub();
}

bool GaussElim::forwardElim() {
    for (int i = 0; i < size; i++){
        int i_max = i;
        int v_max = mat[i_max][i];
        //find the largest pivot
        for (int k = i + 1; k < size; k++){
            if (abs(mat[k][i]) > v_max){
                v_max = mat[k][i], i_max = k;
            }
        }
        //check if the diagonal is zero
        if (!mat[i][i_max]){
            return true; // matrix is singular
        }
        //swap the current row with the max
        if (i_max != i){
            for (int k=0; k<=size; k++){
                double temp = mat[i][k];
                mat[i][k] = mat[i_max][k];
                mat[i_max][k] = temp;
            }
        }
        //matrix manipulation to get the lower triangular to zero
        for (int k= i + 1; k < size; k++){
            double f = mat[k][i] / mat[i][i];
            for (int j= i + 1; j <= size; j++){
                mat[k][j] -= mat[i][j] * f;
            }
            mat[k][i] = 0;
        }
    }
    return false;
}



void GaussElim::backSub() {
    double solution[size];
    for (int i = size - 1; i >= 0; i--){
        solution[i] = mat[i][size];
        for (int j = i+1; j < size; j++){
            solution[i] -= mat[i][j]*solution[j];
        }
        solution[i] = solution[i]/mat[i][i];
    }
    cout << "Solution:" << endl;
    for (int i=0; i<size; i++){
        cout << "x["<<i << "] = " << solution[i] << endl;
    }
}

void GaussElim::print() {
    for(int i = 0; i < size; i++){
        for(int z = 0; z < size+ 1; z++){
            cout << mat[i][z] << " ";
        }
        cout << endl;
    }
}