// additional RHS terms for Cahn-Hillard equation

void add_CH_x(int Nx, int Ny, double dx, double dy, double k, double *C, double *RHS, double *H){

	int i, j, ip1, im1, jp1, jm1;
	double rx, ry;

	rx = k/(dx*dx);
	ry = k/(dy*dy);

	for (i=0; i<Nx; ++i){
		for (j=0; j<Ny; ++j){
			H[i*Ny+j] = pow(C[i*Ny+j], 3) - C[i*Ny+j];
		}
	}

	for (i=0; i<Nx; ++i){

		if (i==0) im1 = Nx-1;
		else im1 = i-1;

		if (i==Nx-1) ip1 = 0;
		else ip1 = i+1;

		for (j=0; j<Ny; ++j){

			if (j==0) jm1 = Ny-1;
			else jm1 = j-1;

			if (j==Ny-1) jp1 = 0;
			else jp1 = j+1;

			// RHS is transposed!
			RHS[j*Nx+i] = RHS[j*Nx+i] + rx*(H[ip1*Ny+j] - 2*H[i*Ny+j] + H[im1*Ny+j])
									  + ry*(H[i*Ny+jp1] - 2*H[i*Ny+j] + H[i*Ny+jm1]);

		}
	}
}


void add_CH_y(int Nx, int Ny, double dx, double dy, double k, double *C, double *RHS, double *H){

	int i, j, ip1, im1, jp1, jm1;
	double rx, ry;

	rx = k/(dx*dx);
	ry = k/(dy*dy);

	for (i=0; i<Nx; ++i){
		for (j=0; j<Ny; ++j){
			H[i*Ny+j] = pow(C[i*Ny+j], 3) - C[i*Ny+j];
		}
	}

	for (i=0; i<Nx; ++i){

		if (i==0) im1 = Nx-1;
		else im1 = i-1;

		if (i==Nx-1) ip1 = 0;
		else ip1 = i+1;

		for (j=0; j<Ny; ++j){

			if (j==0) jm1 = Ny-1;
			else jm1 = j-1;

			if (j==Ny-1) jp1 = 0;
			else jp1 = j+1;

			// no transpose
			RHS[i*Ny+j] = RHS[i*Ny+j] + rx*(H[ip1*Ny+j] - 2*H[i*Ny+j] + H[im1*Ny+j])
									  + ry*(H[i*Ny+jp1] - 2*H[i*Ny+j] + H[i*Ny+jm1]);

		}
	}
}
