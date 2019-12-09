//
// Created by fraktal on 22.11.2019.
//

#ifndef BLUME_CAPEL_METROPOLIS_H
#define BLUME_CAPEL_METROPOLIS_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


mat initial_conf(int L) {
    mat random_conf;
    random_conf.zeros(L, L);

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(-1, 1); // define the range


    for (int i=0; i<L; i++) {
        for(int j=0; j<L; j++) {
            random_conf(i, j) = distr(eng);
        }
    }

    return random_conf;
}

double energy(mat conf, double D, double J) {
    int L = conf.n_rows;
    double part2, part1;
    double value;
    double neig0, neig1, neig2, neig3;
    part1 = 0;
    for (int i=0; i< L; i++) {
        for(int j=0; j<L; j++) {
            value = conf(i, j);
            neig0 = conf((i+1)%L, j);
            neig1 = conf(i, (j+1)%L);
            neig2 = conf((i-1+L)%L, j);
            neig3 = conf(i, (j-1+L)%L);

            part1 += -J*value*(neig0 + neig1 + neig2 + neig3);



        }
    }

    part1 /= 2;
    conf.transform( [](int val) { return (val*val); } );
    part2 = accu(conf);
    part2 *= D;

    return part1 + part2;
}

mat Metropolis(mat initial_conf, double beta, double J, double D, int steps) {
    int L = initial_conf.n_rows;
    mat new_conf;
    new_conf = initial_conf;
    double E0, Enew, delta, r;
    random_device rd;
    mt19937 eng(rd());
    uniform_int_distribution<> distr(0, L-1);
    uniform_int_distribution<> distr2(-1, 1);
    std::uniform_real_distribution<double> distr_double(0.0, 1.0);

    E0 = energy(initial_conf, D, J);
    for(int i=0; i<steps; i++) {
        int site_x = distr(eng);
        int site_y = distr(eng);

        int value = initial_conf(site_x, site_y);
        int new_spin = distr2(eng);

        while(value == new_spin)
            new_spin = distr2(eng);
        new_conf = initial_conf;
        new_conf(site_x, site_y) = new_spin;
        Enew = energy(new_conf, D, J);
        delta = Enew - E0;

        if(delta<0) {
            initial_conf = new_conf;
            E0 = Enew;
        }
        else {
            r = distr_double(eng);
            if(r< exp(-beta*delta)) {
                initial_conf = new_conf;
                E0 = Enew;
            }
        }
    }

    return initial_conf;

}

#endif //BLUME_CAPEL_METROPOLIS_H
