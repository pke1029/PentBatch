
// command: gcc -o hyperdiffusion hyperdiffusion.c -llapack -lblas -O3
// command: ./hyperdiffusion

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "print_matrix.h"
#include "get_penta_mat.h"
#include "get_omega.h"
#include "get_ic.h"
#include "get_RHS.h"

// solve general AX = B
extern void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, 
				   double *B, int *LDB, int *INFO);

// matrix-matrix multiplication
extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, 
				   double *ALPHA, double *A, int *LDA, double *B, 
				   int *LDB, double *BETA, double *C, int *LDC);

// matrix-vector multiplication
extern void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, 
				   double *A, int *LDA, double *X, int *INCX, 
				   double *BETA, double *Y, int *INCY);

int main(){

	int i, j;
	double *temp;

	// parameters **************************************************

	int Nx = 1000;
	double left_boundary = 0.0;
	double right_boundary = 1.0;
	double T = 0.000025;
	double dt = 0.000001;
	
	// *************************************************************

	// parameters for dgemm and dgesv
	char *TRANS, *TRANSA, *TRANSB;
	int M, N, K, LDA, LDB, LDC, NRHS, INCX, INCY, INFO;
	double ALPHA, BETA;
	int *IPIV;
	IPIV = (int*) calloc(Nx, sizeof(int));

	// space grid
	double L = right_boundary - left_boundary;
	double dx = L/Nx;
	double *x_vec;
	x_vec = (double*) calloc(Nx+1, sizeof(double));
	for (i=0; i<Nx+1; ++i) x_vec[i] = left_boundary+i*dx;

	// time grid
	int Nt = floor(T/dt);
	
	// solution storage
	double *sol;
	sol = (double*) calloc((Nt+1)*(Nx+1), sizeof(double));
	
	double *C_old, *dn, *C_new;
	C_old = (double*) calloc(Nx, sizeof(double));
	dn    = (double*) calloc(Nx, sizeof(double));
	C_new = (double*) calloc(Nx, sizeof(double));

	// declare for later use
	double *d_hat, *Ed, *fn;
	d_hat   = (double*) calloc(Nx-2, sizeof(double));
	Ed      = (double*) calloc(Nx-2, sizeof(double));
	fn      = (double*) calloc(4, sizeof(double));

	// diffusivity parameters
	double D = 1.0;
	double gamma = 1.0;
	double r = D*gamma*dt/pow(dx, 4);

	// coefficient
	double a = 0.5*r;
	double b = -2*r;
	double c = 1+3*r;
	double d = -2*r;
	double e = 0.5*r;

	// get pentadiagonal matrix
	N = Nx-2;
	double *E;
	E = (double*) calloc(N*N, sizeof(double));
	get_penta_mat(a, b, c, d, e, N, E);

	// inverse E ***************************************************

	// get idientity matrix
	double *E_inv;
	E_inv = (double*) calloc(N*N, sizeof(double));
	get_penta_mat(0, 0, 1, 0, 0, N, E_inv);

	// parameters for dgesv
	NRHS = N, LDA = N, LDB = N;

	// inverse matrix 
	// warning: E get overwitten by LU factorization on exit, E_inv is E^-1
	dgesv_(&N, &NRHS, E, &LDA, IPIV, E_inv, &LDB, &INFO);

	// *************************************************************

	// check inverse E
	// double *C;
	// C = (double*) calloc(N*N, sizeof(double));
	// get_penta_mat(a, b, c, d, e, N, E);
	// TRANSA = "N", TRANSB = "N", ALPHA = 1, BETA = 0;
	// dgemm_(TRANSA, TRANSB, &N, &N, &N, &ALPHA, E, &N, E_inv, &N, &BETA, C, &N);
	// print_matrix("C", N, N, C);
	// free(C);

	// compute omega = inverse(wee - transpose(g) * inverse(E) * f)
	double *omega;
	omega = (double*) calloc(2*2, sizeof(double));
	get_omega(a, b, c, d, e, Nx-2, E_inv, omega);
	print_matrix("omega", 2, 2, omega);

	// initial condition
	get_ic_cos(Nx, L, x_vec, C_old);
	// get_ic_normal(Nx, L, x_vec, C_old);

	// write to solution
	for (i=0; i<Nx; ++i) sol[i] = C_old[i];
	sol[Nx] = C_old[0];
	
	for (j=0; j<Nt; ++j){

		get_RHS(Nx, r, C_old, dn);
		for (i=0; i<Nx-2; ++i) d_hat[i] = dn[i];

		// compute Ed = inverse(E) * d_hat
		N = Nx-2;
		TRANS = "N", M = N, N = N, ALPHA = 1.0, LDA = N, INCX = 1, BETA = 0, INCY = 1;
		dgemv_(TRANS, &M, &N, &ALPHA, E_inv, &LDA, d_hat, &INCX, &BETA, Ed, &INCY);

		// lambda = small_d * - transpose(g) * inverse(E) *d_hat
		dn[Nx-2] = dn[Nx-2] - (e*Ed[0] + a*Ed[N-2] + b*Ed[N-1]);
		dn[Nx-1] = dn[Nx-1] - (d*Ed[0] + e*Ed[1]   + a*Ed[N-1]);

		// small_x = omega * lambda
		C_new[Nx-2] = omega[0]*dn[Nx-2] + omega[1]*dn[Nx-1];
		C_new[Nx-1] = omega[2]*dn[Nx-2] + omega[3]*dn[Nx-1];

		// compute non-zero entries of f * small_x
		fn[0] = a*C_new[Nx-2] + b*C_new[Nx-1];
		fn[1] = a*C_new[Nx-1];
		fn[2] = e*C_new[Nx-2];
		fn[3] = d*C_new[Nx-2] + e*C_new[Nx-1];

		// X_hat = inverse(E) * d_hat - inverse(E) * f * small_x
		for (i=0; i<N; ++i){

			C_new[i] = Ed[i] - (fn[0]*E_inv[i*N+0] 
								+ fn[1]*E_inv[i*N+1] 
								+ fn[2]*E_inv[i*N+(N-2)] 
								+ fn[3]*E_inv[i*N+(N-1)]);
		}

		// write to solution
		for (i=0; i<Nx; ++i) sol[(j+1)*(Nx+1)+i] = C_new[i];
		sol[(j+1)*(Nx+1)+Nx] = C_new[0];

		// update
		temp = C_new;
        C_new = C_old;
        C_old = temp;
	}

	// write solution to file **************************************

    FILE *fptr;
    fptr = fopen("plot_hyperdiffusion.txt", "w");
    if(fptr == NULL){
        printf("error opening file");
        exit(1);
    }

    for (j=0; j<(Nt+1); ++j){
    	for (i=0; i<(Nx+1); ++i){

    		if (i!=Nx) fprintf(fptr, "%f,", sol[j*(Nx+1)+i]);
            else fprintf(fptr, "%f\n", sol[j*(Nx+1)+i]);
    	}
    }

    fclose(fptr);

    // *************************************************************
	
	free(IPIV);
	free(x_vec);
	free(sol);
	free(C_old);
	free(dn);
	free(C_new);
	free(d_hat);
	free(Ed);
	free(fn);
	free(E);
	free(E_inv);
	free(omega);

	return 0;
}