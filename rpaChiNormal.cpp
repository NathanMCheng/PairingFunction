/*
** main.c
*/
#include <iostream>
#include <complex>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

typedef complex<double> dcomp;

const dcomp I(0.0, 1.0);
const double pi = 3.1415926535897;

int main (int argc, char* argv[]) { //input frequency size, qx size, qy size

    int wsize=0;
    double wstart = 0.;
    double Gamma = 0.;
    int meshSize = 0;
    string name;

    ofstream outputchi;



    if(argc == 6) {
        wsize = atoi(argv[1]);
        wstart = atof(argv[2]);
        Gamma = atof(argv[3]);
        meshSize = atoi(argv[4]);
        name = argv[5];
    }

    outputchi.open(name.c_str());

    //    dcomp **chi;
    //    chi = (dcomp **) malloc(matrixsize1*sizeof(dcomp*));
    //    for(int i = 0; i < matrixsize1; i++) {
    //        chi[i] = (dcomp *) malloc(matrixsize2*sizeof(dcomp));
    //    }

    double delta0,chemPotent,t1,t2,t3,t4,t5,U,KpereV,beta,N, qxi, qyi, kxi, kyi, cosKx, cosKy, cos2Kx, cos2Ky,
            cosKQx, cosKQy, cos2KQx, cos2KQy, epsilonK, deltaK, epsilonKQ, deltaKQ, EK, EKQ, CplusQK, CminusQK,
            fermiEKQ, fermiEK,chi,w; //fermiMinusEKQ, fermiMinusEK,
    dcomp sum, chi0;

    // Table Variables
    delta0 = 0.035; //eV (35.0 meV) superconducting gap
    chemPotent = 0.0879; //eV (87.9 meV)
    t1 = -0.5547; //eV (-554.7 meV)
    t2 = 0.1327; //eV (132.7 meV)
    t3 = 0.0132; //eV (13.2 meV)
    t4 = -0.1849; //eV (184.9 meV)
    t5 = 0.0265; //eV (26.5 meV)
    U = 0.165; //%eV (165.0 meV) U-limit of the Hubbard model

    // %additional variables needed
    KpereV = 11604.505;
    beta = 1./(300./KpereV); /*%1/eV*/
    //Gamma = 0.0024; /*%eV (1 meV)*/

    N = ((double) meshSize)*((double) meshSize); /*% number of sites on a square lattice wiht periodic boundary conditions*/
    for (int z = 0; z<wsize; z++) {
        w = (wstart+ ( (double) z))/1000.;
        for (int i = 0; i < meshSize+1; i++) {
            qxi = (0.5+( (double) i)/( (double) meshSize))*pi;
            for (int j = 0; j < meshSize+1; j++) {
                qyi = (0.5+( (double) j)/( (double) meshSize))*pi;
                sum = 0.+0.*I;
                for (int k = 0; k <meshSize+1; k++) {
                    kxi = (double (-1.+2.*( (double) k)/( (double) meshSize)))*pi;
                    for (int l = 0; l<meshSize+1; l++) {
                        kyi = (double (-1.+2.*( (double) l)/( (double) meshSize)))*pi;

                        cosKx = cos(kxi);
                        cosKy = cos(kyi);
                        cos2Kx = cos(2.*kxi);
                        cos2Ky = cos(2.*kyi);

                        cosKQx = cos(qxi+kxi);
                        cosKQy = cos(qyi+kyi);
                        cos2KQx = cos(2.*(qxi+kxi));
                        cos2KQy = cos(2.*(qyi+kyi));

                        epsilonK =  t1/2.*(cosKx+cosKy)+t2*cosKx*cosKy +t3/2.*(cos2Kx+cos2Ky) +
                                t4/2.*(cos2Kx*cosKy+cosKx*cos2Ky)+t5*cos2Kx*cos2Ky+chemPotent;
                        //band parameters at various points in the BZ

                        // deltaK = delta0/2.*(cosKx-cosKy); //gap at various points in the BZ

                        epsilonKQ =  t1/2.*(cosKQx+cosKQy)+t2*cosKQx*cosKQy +t3/2.*(cos2KQx+cos2KQy) +
                                t4/2.*(cos2KQx*cosKQy+cosKQx*cos2KQy)+t5*cos2KQx*cos2KQy+chemPotent;
                        //band parameters at various points in the BZ

                        //deltaKQ = delta0/2.*(cosKQx-cosKQy); //gap at various points in the BZ

                        //EK = sqrt(pow(epsilonK,2)+pow(deltaK,2)); //BCS dispersion of quasiparticles

                        //EKQ = sqrt(pow(epsilonKQ,2)+pow(deltaKQ,2)); //BCS dispersion of qusiparticles

                        //CplusQK = (1. + (epsilonKQ*epsilonK+deltaKQ*deltaK)/(EKQ*EK)); //coherence factor

                        //CminusQK = (1. - (epsilonKQ*epsilonK+deltaKQ*deltaK)/(EKQ*EK)); //coherence factor

                        fermiEKQ = 1./(exp(beta*epsilonKQ)+1.);

                        fermiEK = 1./(exp(beta*epsilonK)+1.);

                        //                        fermiMinusEK = 1./(exp(-beta*EK)+1.);

                        //                        fermiMinusEKQ = 1./(exp(-beta*EKQ)+1.);

                        sum = sum + (fermiEKQ-fermiEK)/(w+Gamma*I-(epsilonKQ-epsilonK));
                        //                        sum =   sum - CplusQK*(fermiEKQ-fermiEK)/(w+Gamma*I+EKQ-EK) +
                        //                                CminusQK*(fermiEKQ-fermiMinusEK)/(w+Gamma*I-EKQ-EK) +
                        //                                CminusQK*(fermiMinusEKQ-fermiEK)/(w+Gamma*I+EKQ+EK) -
                        //                                CplusQK*(fermiMinusEKQ+fermiMinusEK)/(w+Gamma*I-EKQ+EK);


                    }

                }
                chi0 = sum/N;
                chi = imag(chi0)/(pow(1.-U*(real(chi0)),2)+pow(U,2)*pow(imag(chi0),2));
                outputchi << chi << "\t";
            }
            outputchi << "\n";
        }
        //outputchi <<"\r";
    }

    outputchi.close();
    ofstream finishfile;
    string finishname = "job_finished";

    finishfile.open(finishname.c_str());
    finishfile.close();
}

