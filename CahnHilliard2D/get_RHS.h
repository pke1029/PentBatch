
void get_RHS_X(int Nx, int Ny, double ry, double rxy, double *C_old, double *RHS){

	int i, j;

	int ip1, im1, jp2, jp1, jm1, jm2;

	for (i=0; i<Nx; ++i){

		if (i==0) im1 = Nx-1;
		else im1 = i-1;

		if (i==Nx-1) ip1 = 0;
		else ip1 = i+1;

		for(j=0; j<Ny; ++j){

			if (j==0) jm1 = Ny-1, jm2 = Ny-2;
			else if (j==1) jm1 = 0, jm2 = Ny-1;
			else jm1 = j-1 , jm2 = j-2;

			if (j==Ny-1) jp1 = 0, jp2 = 1;
			else if (j==Ny-2) jp1 = Ny-1, jp2 = 0;
			else jp1 = j+1, jp2 = j+2;

			// RHS is transposed!!!

			// 3 by 3 cross term
			RHS[j*Nx+i] = rxy * (    C_old[im1*Ny+jm1] - 2*C_old[im1*Ny+j] +   C_old[im1*Ny+jp1]
								 - 2*C_old[i*Ny+jm1]   + 4*C_old[i*Ny+j]   - 2*C_old[i*Ny+jp1]
								 +   C_old[ip1*Ny+jm1] - 2*C_old[ip1*Ny+j] +   C_old[ip1*Ny+jp1]);

			// (add onto it) 1 by 5 dyyyy term
			RHS[j*Nx+i] = C_old[i*Ny+j] - RHS[j*Nx+i] 
						  - 0.5*ry*(C_old[i*Ny+jm2] - 4*C_old[i*Ny+jm1] + 6*C_old[i*Ny+j]
									- 4*C_old[i*Ny+jp1] + C_old[i*Ny+jp2]);
		}
	}	
}


void get_RHS_Y(int Nx, int Ny, double rx, double rxy, double *C_half, double *RHS){

	int i, j;

	int ip2, ip1, im1, im2, jp1, jm1;

	for (i=0; i<Nx; ++i){

		if (i==0) im1 = Nx-1, im2 = Nx-2;
		else if (i==1) im1 = 0, im2 = Nx-1;
		else im1 = i-1, im2 = i-2;

		if (i==Nx-1) ip1 = 0, ip2 = 1;
		else if (i==Nx-2) ip1 = Nx-1, ip2 = 0;
		else ip1 = i+1, ip2 = i+2;

		for(j=0; j<Ny; ++j){

			if (j==0) jm1 = Ny-1;
			else jm1 = j-1;

			if (j==Ny-1) jp1 = 0;
			else jp1 = j+1;

			// 3 by 3 cross term
			RHS[i*Ny+j] = rxy * (    C_half[im1*Ny+jm1] - 2*C_half[im1*Ny+j] +   C_half[im1*Ny+jp1]
								 - 2*C_half[i*Ny+jm1]   + 4*C_half[i*Ny+j]   - 2*C_half[i*Ny+jp1]
								 +   C_half[ip1*Ny+jm1] - 2*C_half[ip1*Ny+j] +   C_half[ip1*Ny+jp1]);

			// (add onto it) 1 by 5 dxxxx term
			RHS[i*Ny+j] = C_half[i*Ny+j] - RHS[i*Ny+j] 
						  - 0.5*rx*(C_half[im2*Ny+j] - 4*C_half[im1*Ny+j] + 6*C_half[i*Ny+j]
									- 4*C_half[ip1*Ny+j] + C_half[ip2*Ny+j]);
		}
	}	
}


void get_RHS_stable(int Nx, int Ny, double rx, double ry, double rxy, double *C_old, 
					double *C_new, double *C_bar, double *RHS){

	int i, j;
	double dxxxx, dxxyy, dyyyy;

	int ip2, ip1, im1, im2, jp2, jp1, jm1, jm2;

	for (i=0; i<Nx; ++i){

		if (i==0) im1 = Nx-1, im2 = Nx-2;
		else if (i==1) im1 = 0, im2 = Nx-1;
		else im1 = i-1, im2 = i-2;

		if (i==Nx-1) ip1 = 0, ip2 = 1;
		else if (i==Nx-2) ip1 = Nx-1, ip2 = 0;
		else ip1 = i+1, ip2 = i+2;

		for(j=0; j<Ny; ++j){

			if (j==0) jm1 = Ny-1, jm2 = Ny-2;
			else if (j==1) jm1 = 0, jm2 = Ny-1;
			else jm1 = j-1 , jm2 = j-2;

			if (j==Ny-1) jp1 = 0, jp2 = 1;
			else if (j==Ny-2) jp1 = Ny-1, jp2 = 0;
			else jp1 = j+1, jp2 = j+2;

			// 3 by 3 dxxyy term
			dxxyy = 2*rxy * (    C_bar[im1*Ny+jm1] - 2*C_bar[im1*Ny+j] +   C_bar[im1*Ny+jp1]
						     - 2*C_bar[i*Ny+jm1]   + 4*C_bar[i*Ny+j]   - 2*C_bar[i*Ny+jp1]
							 +   C_bar[ip1*Ny+jm1] - 2*C_bar[ip1*Ny+j] +   C_bar[ip1*Ny+jp1]);

			// dyyyy term
			dyyyy = ry * (C_bar[i*Ny+jm2] - 4*C_bar[i*Ny+jm1] + 6*C_bar[i*Ny+j]
							- 4*C_bar[i*Ny+jp1] + C_bar[i*Ny+jp2]);

			// dxxxx term
			dxxxx = rx * (C_bar[im2*Ny+j] - 4*C_bar[im1*Ny+j] + 6*C_bar[i*Ny+j]
							- 4*C_bar[ip1*Ny+j] + C_bar[ip2*Ny+j]);

			// RHS is transposed for memory access!
			RHS[j*Nx+i] = -2.0/3.0 * (C_new[i*Ny+j] - C_old[i*Ny+j]) - 2.0/3.0 * (dxxxx + dxxyy + dyyyy);
		}
	}
}
