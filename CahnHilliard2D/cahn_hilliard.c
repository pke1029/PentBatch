
// command: gcc -o cahn_hilliard cahn_hilliard.c -lm -O3
// command: ./cahn_hilliard

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "print_matrix.h"
#include "matrix_multiplication.h"
#include "get_ic.h"
#include "get_LHS.h"
#include "solve_penta.h"
#include "get_omega.h"
#include "get_RHS.h"
#include "solve_cyclic_penta.h"
#include "add_CH.h"


int main(){

	int i, j, t;

	// parameters **************************************************
	
	int Nx = 64;
	int Ny = 64;
	double T = 2.0;
	double dt = 0.0001;
	double dt_new = 0.0001;
	double left_boundary = 0.0;
	double right_boundary = 2.0 * M_PI;
	double bottom_boundary = 0.0;
	double top_boundary = 2.0 * M_PI;

	// *************************************************************
	
	double *tempPtr;

	// grid
	double Lx = right_boundary - left_boundary;
	double Ly = top_boundary - bottom_boundary;
	double dx = Lx/Nx;
	double dy = Ly/Ny;
	int Nt = floor((T-dt)/dt_new) + 1;
	double *x_vec, *y_vec;
	x_vec = (double*) calloc(Nx, sizeof(double));
	y_vec = (double*) calloc(Ny, sizeof(double));
	for (i=0; i<Nx; ++i) x_vec[i] = left_boundary + i*dx;
	for (i=0; i<Ny; ++i) y_vec[i] = bottom_boundary + i*dy;

	// old and new matrix
	double *C_old, *C_half, *C_new, *C_bar, *H;
	C_old = (double*) calloc(Nx*Ny, sizeof(double));
	C_half = (double*) calloc(Nx*Ny, sizeof(double));
	C_new = (double*) calloc(Nx*Ny, sizeof(double));
	C_bar = (double*) calloc(Nx*Ny, sizeof(double));
	H = (double*) calloc(Nx*Ny, sizeof(double));

	// declare for later use
	double *RHS, *dn_x, *dn_y;
	RHS = (double*) calloc(Nx*Ny, sizeof(double));

	// coefficient
	double D = 0.1;
	double gamma = 0.01;
	double rx = D*gamma*dt/pow(dx, 4);
	double ry = D*gamma*dt/pow(dy, 4);
	double rxy = D*gamma*dt/(dx*dx*dy*dy);
	double k = D*dt/2.0;
	double a, b, c, d, e;
	
	// initialize matrix
	double *alpha_x, *beta_x, *epsilon_x, *gamma_x, *delta_x;
	double *alpha_y, *beta_y, *epsilon_y, *gamma_y, *delta_y;
	alpha_x   = (double*) calloc(Nx, sizeof(double));
	beta_x    = (double*) calloc(Nx, sizeof(double));
	epsilon_x = (double*) calloc(Nx, sizeof(double));
	gamma_x   = (double*) calloc(Nx, sizeof(double));
	delta_x   = (double*) calloc(Nx, sizeof(double));
	alpha_y   = (double*) calloc(Ny, sizeof(double));
	beta_y    = (double*) calloc(Ny, sizeof(double));
	epsilon_y = (double*) calloc(Ny, sizeof(double));
	gamma_y   = (double*) calloc(Ny, sizeof(double));
	delta_y   = (double*) calloc(Ny, sizeof(double));

	double *f_x, *g_x, *omega_x;
	double *f_y, *g_y, *omega_y;
	f_x  = (double*) calloc((Nx-2)*2, sizeof(double));
	g_x = (double*) calloc(2*(Nx-2), sizeof(double));
	omega_x = (double*) calloc(4, sizeof(double));
	f_y  = (double*) calloc((Ny-2)*2, sizeof(double));
	g_y = (double*) calloc(2*(Ny-2), sizeof(double));
	omega_y = (double*) calloc(4, sizeof(double));

	// temporary storage for solve_cyclic_penta
	double *small_d;
	small_d = (double*) calloc(2, sizeof(double));

	// get initial condition
	// get_sin_ic(Nx, Ny, Lx, Ly, x_vec, y_vec, C_old);
	get_rand_ic(Nx, Ny, Lx, Ly, x_vec, y_vec, C_old);

	// *************************************************************
	// conditionally stable scheme                                 *
	// *************************************************************

	// start clock
	clock_t start = clock(), diff;

	// omega_x *****************************************************
	
	// get E_x, f_x, gT_x
	int Nx2 = Nx - 2;
	a = 0.5*rx, b = -2*rx, c = 1 + 3*rx, d = -2*rx, e = 0.5*rx;
	get_penta_mat(Nx2, a, b, c, d, e, epsilon_x, beta_x, alpha_x, 
				  gamma_x, delta_x);
	get_fg(Nx2, a, b, d, e, f_x, g_x);

	// LU factorize E_x
	LU_factorisation(Nx2, epsilon_x, beta_x, alpha_x, gamma_x, 
					 delta_x);
	// compute transpose(g)*inverse(E) {overwrite g on out}
	solve_pentaX(Nx2, 2, alpha_x, beta_x, epsilon_x, gamma_x, 
				 delta_x, g_x);
	// compute omega
	get_omega(Nx2, b, c, d, alpha_x, beta_x, epsilon_x, gamma_x, 
			  delta_x, g_x, f_x, omega_x);
	
	// *************************************************************

	// omega_y *****************************************************
	
	// get E_y, f_y, gT_y
	int Ny2 = Ny - 2;
	a = 0.5*ry, b = -2*ry, c = 1 + 3*ry, d = -2*ry, e = 0.5*ry;
	get_penta_mat(Ny2, a, b, c, d, e, epsilon_y, beta_y, alpha_y, 
				  gamma_y, delta_y);
	get_fg(Ny2, a, b, d, e, f_y, g_y);

	// LU factorize E_y
	LU_factorisation(Ny2, epsilon_y, beta_y, alpha_y, gamma_y, 
					 delta_y);
	// compute transpose(g)*inverse(E) {overwrite g on out}
	solve_pentaX(Ny2, 2, alpha_y, beta_y, epsilon_y, gamma_y, 
				 delta_y, g_y);
	// compute omega
	get_omega(Ny2, b, c, d, alpha_y, beta_y, epsilon_y, gamma_y, 
			  delta_y, g_y, f_y, omega_y);

	// *************************************************************

	// iterate for 1 timestep **************************************

	get_RHS_X(Nx, Ny, ry, rxy, C_old, RHS);
	add_CH_x(Nx, Ny, dx, dy, k, C_old, RHS, H);

	for (j=0; j<Ny; ++j){

		dn_x = &RHS[j*Nx];
		
		solve_cyclic_penta(Nx, alpha_x, beta_x, epsilon_x, gamma_x, 
						   delta_x, f_x, g_x, omega_x, dn_x, 
						   small_d);

		// undo transpose
		for (i=0; i<Nx; ++i) C_half[i*Ny+j] = dn_x[i];
	}

	get_RHS_Y(Nx, Ny, rx, rxy, C_half, RHS);
	add_CH_y(Nx, Ny, dx, dy, k, C_half, RHS, H);

	for (i=0; i<Nx; ++i){

		dn_y = &RHS[i*Ny];

		solve_cyclic_penta(Ny, alpha_y, beta_y, epsilon_y, gamma_y, 
						   delta_y, f_y, g_y, omega_y, dn_y, 
						   small_d);
	}

	tempPtr = C_new;
	C_new = RHS;
	RHS = tempPtr;

	// *************************************************************
	
	// *************************************************************
	// unconditionally stable scheme                               *
	// *************************************************************

	dt = dt_new;
	rx = D*gamma*dt/pow(dx, 4);
	ry = D*gamma*dt/pow(dy, 4);
	rxy = D*gamma*dt/(dx*dx*dy*dy);
	k = 2.0/3.0 * dt * D;

	// omega_x *****************************************************

	// get E_x, f_x, gT_x
	a = 2*rx/3, b = -8*rx/3, c = 1 + 4*rx, d = -8*rx/3, e = 2*rx/3;
	get_penta_mat(Nx2, a, b, c, d, e, epsilon_x, beta_x, alpha_x, 
				  gamma_x, delta_x);
	get_fg(Nx2, a, b, d, e, f_x, g_x);

	// LU factorize E_x
	LU_factorisation(Nx2, epsilon_x, beta_x, alpha_x, gamma_x, 
					 delta_x);
	// compute transpose(g)*inverse(E) {overwrite g on out}
	solve_pentaX(Nx2, 2, alpha_x, beta_x, epsilon_x, gamma_x, 
				 delta_x, g_x);
	// compute omega
	get_omega(Nx2, b, c, d, alpha_x, beta_x, epsilon_x, gamma_x, 
			  delta_x, g_x, f_x, omega_x);

	// *************************************************************

	// omega_y *****************************************************

	// get E_y, f_y, gT_y
	a = 2*ry/3, b = -8*ry/3, c = 1 + 4*ry, d = -8*ry/3, e = 2*ry/3;
	get_penta_mat(Ny2, a, b, c, d, e, epsilon_y, beta_y, alpha_y, 
				  gamma_y, delta_y);
	get_fg(Ny2, a, b, d, e, f_y, g_y);

	// LU factorize E_y
	LU_factorisation(Ny2, epsilon_y, beta_y, alpha_y, gamma_y, 
					 delta_y);
	// compute transpose(g)*inverse(E) {overwrite g on out}
	solve_pentaX(Ny2, 2, alpha_y, beta_y, epsilon_y, gamma_y, 
				 delta_y, g_y);
	// compute omega
	get_omega(Ny2, b, c, d, alpha_y, beta_y, epsilon_y, gamma_y, 
			  delta_y, g_y, f_y, omega_y);

	// *************************************************************
	
	for (t=1; t<Nt; ++t){

		for (i=0; i<Nx*Ny; ++i) C_bar[i] = 2*C_new[i] - C_old[i];

		get_RHS_stable(Nx, Ny, rx, ry, rxy, C_old, C_new, C_bar, 
					   RHS);
		add_CH_x(Nx, Ny, dx, dy, k, C_new, RHS, H);

		for (j=0; j<Ny; ++j){

			dn_x = &RHS[j*Nx];

			solve_cyclic_penta(Nx, alpha_x, beta_x, epsilon_x, 
							   gamma_x, delta_x, f_x, g_x, omega_x, 
							   dn_x, small_d);

			// undo transpose
			for (i=0; i<Nx; ++i) C_half[i*Ny+j] = dn_x[i];
		}

		for (i=0; i<Nx; ++i){

			dn_y = &C_half[i*Ny];

			solve_cyclic_penta(Ny, alpha_y, beta_y, epsilon_y, 
							   gamma_y, delta_y, f_y, g_y, omega_y, 
							   dn_y, small_d);
		}

		// overwrite C_old for C_newest
		for (i=0; i<Nx*Ny; ++i) C_old[i] = C_half[i] + C_bar[i];

		// swap C_old with C_new
		tempPtr = C_old;
		C_old = C_new;
		C_new = tempPtr;
	}

	// stop timer
	diff = clock() - start;
	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds", msec/1000, 
		   msec%1000);
	printf("\n");

	// write to file
	FILE *fptr;
    fptr = fopen("plot_cahn_hilliard.txt", "w");
    for (i=0; i<Nx; ++i){
    	for (j=0; j<Ny; ++j){
    		if (j!=Ny-1) fprintf(fptr, "%f,", C_new[i*Ny+j]);
            else fprintf(fptr, "%f\n", C_new[i*Ny+j]);
    	}
    }
    fclose(fptr);

    free(small_d);
    free(x_vec);
    free(y_vec);
    free(C_old);
    free(C_half);
    free(C_new);
    free(C_bar);
    free(H);
    free(RHS);
    free(alpha_x);
    free(alpha_y);
    free(beta_x);
    free(beta_y);
    free(epsilon_x);
    free(epsilon_y);
    free(gamma_x);
    free(gamma_y);
    free(delta_x);
    free(delta_y);
    free(f_x);
    free(f_y);
    free(g_x);
    free(g_y);
    free(omega_x);
    free(omega_y);

	return 0;
}
