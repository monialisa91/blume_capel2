#include <iostream>
#include <armadillo>
#include "metropolis.h"
#include "observables.h"
#include "wolff.h"

using namespace std;
using namespace arma;

// BLUME-CAPEL MODEL

int main() {

    int L = 16;
    mat initial = initial_conf(L);
    double J = 1;
    double D = 0.0;
    double T, beta;
    int mes = 10000;
    vec energies(mes);
    vec magn(mes);
    double E, M;

    mat configurations;
    mat new_conf;
//    configurations.load("L=10_T=1.51.txt");
//    new_conf = configurations.row(configurations.n_rows - 1);
//    new_conf.reshape(L, L);
//    initial = new_conf;
    initial = initial_conf(L);

    ostringstream ss1;
    ss1 << L;
    string Size = ss1.str();

    ostringstream ssD;
    ssD << D;
    string D_par = ssD.str();

    string conf_name, hist_energy_name, hist_magn_name, capacity_name, suscpet_name;

    suscpet_name = "L=" + Size + "_suscp.txt";
    capacity_name = "L=" + Size + "_heat.txt";


    //1. Thermalisation

    for(int i=230; i>130; i--) {
        T = i*0.01;
        cout << T << endl;
        beta = 1.0/T;

        ostringstream ss2;

        ss2 << T;
        string Temp = ss2.str();

        conf_name = "L=" + Size  + "D=" + D_par + "_T=" + Temp + ".txt";
        hist_energy_name = "L=" + Size  + "D=" + D_par +  "_T=" + Temp + "_energy.txt";
        hist_magn_name = "L=" + Size + "D=" + D_par +"_T=" + Temp + "_magnet.txt";

        for(int m=0; m<2000; m++) {
            initial = Metropolis(initial, beta, J, D, 3*L);
            initial = WolffUpdate(initial, beta, J, 1);
        }
        energies.zeros();
        magn.zeros();

        ofstream myfile;
        myfile.open(conf_name.c_str(), ios_base::app);

        //2. MEASUREMENT
        for(int k=0; k<mes; k++) {
            for(int m=0; m<L; m++) {
                initial = Metropolis(initial, beta, J, D, 3*L);
                initial = WolffUpdate(initial, beta, J, 1);
            }
            E = energy(initial, D, J);
            M  = magne(initial);
            energies(k) = E;
            magn(k) = M;
            saveintofile(initial, myfile);
        }
        myfile.close();
        cout << "heat_capa " << heat_capacity(beta, energies) << endl;
        energies.save(hist_energy_name, csv_ascii);
        magn.save(hist_magn_name, csv_ascii);

    }

    return 0;
}