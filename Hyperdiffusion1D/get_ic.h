
#include <math.h>


void get_ic_cos(int Nx, double L, double *x_vec, double *fn){

	int i;

	int n = 2;
	double amplitude = 1.0;

	for (i=0; i<Nx; ++i){
		fn[i] = amplitude*cos(n*2*M_PI*x_vec[i]/L);
	}
}


void get_ic_normal(int Nx, double L, double *x_vec, double *fn){

	int i;

	double amplitude = 1.0;

	for (i=0; i<Nx; ++i){
        fn[i] = amplitude*exp(-20*pow((x_vec[i]-L/2),2));
    }
}


void get_ic_step(int Nx, double L, double *x_vec, double *fn){

	int = i;

	double amplitude = 1.0;

	for (i=0; i<Nx; ++i){
		if (i<Nx/2) fn[i] = amplitude;
        else fn[i] = 0;
    }
}