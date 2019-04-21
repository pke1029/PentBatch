
void get_sin_ic(int Nx, int Ny, double Lx, double Ly, double *x_vec, double *y_vec, double *C_old){

	int i, j;
	double temp1, temp2;

	int n = 1;
	double amplitude = 1.0;

	for (i=0; i<Nx; ++i){

		temp1 = sin(n*2*M_PI*x_vec[i]/Lx);

		for (j=0; j<Ny; ++j){

			temp2 = sin(n*2*M_PI*y_vec[j]/Ly);

			C_old[i*Ny+j] = amplitude * temp1 * temp2;
		}
	}
}

void get_rand_ic(int Nx, int Ny, double Lx, double Ly, double *x_vec, double *y_vec, double *C_old){

	int i, j;

	double a = -0.1;
	double b = 0.1;

	for (i=0; i<Nx; ++i){
		for (j=0; j<Ny; ++j){
			C_old[i*Ny+j] = (b - a) * (double)rand() / (double)RAND_MAX + a;
		}
	}
}
