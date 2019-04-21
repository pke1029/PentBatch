
void get_RHS(int Nx, double r, double *fn, double *dn){

	int i;

	double a, b, c, d, e;
	a = -0.5*r;
	b = 2*r;
	c = 1-3*r;
	d = 2*r;
	e = -0.5*r;

	dn[0] = a*fn[Nx-2] + b*fn[Nx-1] + c*fn[0] + d*fn[1] + e*fn[2];
	dn[1] = a*fn[Nx-1] + b*fn[0] + c*fn[1] + d*fn[2] + e*fn[3];

	for (i=2; i<Nx-2; ++i){

		dn[i] = a*fn[i-2] + b*fn[i-1] + c*fn[i] + d*fn[i+1] + e*fn[i+2];
	}

	dn[Nx-2] = a*fn[Nx-4] + b*fn[Nx-3] + c*fn[Nx-2] + d*fn[Nx-1] + e*fn[0];
	dn[Nx-1] = a*fn[Nx-3] + b*fn[Nx-2] + c*fn[Nx-1] + d*fn[0] + e*fn[1];
}
