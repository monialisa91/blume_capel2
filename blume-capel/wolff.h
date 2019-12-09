//
// Created by fraktal on 27.11.2019.
//

#ifndef BLUME_CAPEL_WOLFF_H
#define BLUME_CAPEL_WOLFF_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat WolffUpdate(mat conf, double beta, double J, int N) {
    int L = conf.n_rows;
    random_device rd;
    mt19937 eng(rd());
    uniform_int_distribution<> distr(0, L-1);

    vec sites_x(L*L);
    vec sites_y(L*L);

    vec neig_sites_x(4);
    vec neig_sites_y(4);

    cube neigb (4, L, L);
    mat already_cluster(L, L);


    int ss;
    int x_site, y_site;
    double r, p;





    // Number of iterations
    for(int it=0; it<N; it++) {

        for(int i=0; i<L; i++) {
            for(int j=0; j<L; j++) {
                neigb(0, i, j) = conf(i, (j+1)%L);
                neigb(1, i, j) = conf((i+L-1)%L, j);
                neigb(2, i, j) = conf(i, (j+L-1)%L);
                neigb(3, i, j) = conf((i+1)%L, j);
            }
        }



//        cout << it << endl;

        sites_x.zeros();
        sites_y.zeros();

        sites_x(0) = distr(eng);
        sites_y(0) = distr(eng);


        ss = conf(sites_x(0), sites_y(0));

        while(ss==0) {
            sites_x(0) = distr(eng);
            sites_y(0) = distr(eng);
            ss = conf(sites_x(0), sites_y(0));
        }


        already_cluster.zeros();
        already_cluster(sites_x(0), sites_y(0)) = 1;


        int k = 1; // number of sites in the cluster
        int limit0 = 0;
        int limit1 = k;

        while (limit0 < limit1) {

            for (int i = limit0; i < limit1; i++) {
                x_site = (int) sites_x(i);
                y_site = (int) sites_y(i);

                neig_sites_x(0) = x_site;
                neig_sites_x(1) = (x_site + L - 1) % L;
                neig_sites_x(2) = x_site;
                neig_sites_x(3) = (x_site + 1) % L;

                neig_sites_y(0) = (y_site + 1) % L;
                neig_sites_y(1) = y_site;
                neig_sites_y(2) = (y_site - 1 + L) % L;
                neig_sites_y(3) = y_site;



                for (int j = 0; j < 4; j++) {
                    if (neigb(j, sites_x(i), sites_y(i)) == ss &&
                        already_cluster(neig_sites_x(j), neig_sites_y(j)) != 1.0) {
                        r = ((double) rand() / (RAND_MAX));
                        p = 1-exp(-2 * beta);
                        if (r < p) {
                            sites_x(k) = neig_sites_x(j);
                            sites_y(k) = neig_sites_y(j);
                            already_cluster(neig_sites_x(j), neig_sites_y(j)) = 1;
                            k++;

                        }

                    }
                }

            }


            limit0 = limit1;
            limit1 = k;

        }





//        sites_x = sites_x.head(k);
//        sites_y = sites_y.head(k);


        int new_s = -1*ss;



        for (int n = 0; n < k; n++) {
            conf(sites_x(n), sites_y(n)) = new_s;
        }

//        conf.print();

    }
    return conf;

}

#endif //BLUME_CAPEL_WOLFF_H
