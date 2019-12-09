//
// Created by fraktal on 27.11.2019.
//

#ifndef BLUME_CAPEL_OBSERVABLES_H
#define BLUME_CAPEL_OBSERVABLES_H

#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;

double magne(mat conf) {
    int L = conf.n_rows;
    double suma = 0;
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++) {
            suma += conf(i, j);
        }
    }

    return suma;
}

double heat_capacity(double beta, vec energy) {
    int size = energy.n_rows;
    double suma = accu(energy)/size;
    double suma_kw = accu(energy%energy)/size;

    return beta*beta*(suma_kw - suma*suma);
}

double susceptibility(double beta, vec magn) {
    int size = magn.n_rows;
    double suma = accu(magn);
    double suma_kw = accu(magn%magn);

    return beta*(suma_kw - suma*suma);
}

void saveintofile(mat config, ofstream &file) {
    int L = config.n_rows;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            file << setprecision(1) << config(i, j) << " ";
        }
    }
    file << "\n";
}


#endif //BLUME_CAPEL_OBSERVABLES_H
