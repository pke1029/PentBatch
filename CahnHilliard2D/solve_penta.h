
// solve Ax = b

void LU_factorisation(int N, double *a, double *b, double *c, double *d, double *e){

	// decompose A into L*U
	
	// on out: a_vec become epsilon
	//		   b_vec become beta
	//		   c_vec become alpha
	//		   d_vec become gamma
	// 		   e_vec become delta

	int i;

	c[0] = c[0];
	d[0] = d[0]/c[0];
	e[0] = e[0]/c[0];
	
	b[1]  = b[1];
	c[1] = c[1]-b[1]*d[0];
	d[1] = (d[1]-b[1]*e[0])/c[1];
	e[1] = e[1]/c[1];

	for (i=2; i<N-2; ++i){
		b[i]  = b[i]-a[i]*d[i-2];
		c[i] = c[i]-a[i]*e[i-2]-b[i]*d[i-1];
		d[i] = (d[i]-b[i]*e[i-1])/c[i];
		e[i] = e[i]/c[i];
	}

	b[N-2]  = b[N-2]-a[N-2]*d[N-4];
	c[N-2] = c[N-2]-a[N-2]*e[N-4]-b[N-2]*d[N-3];
	d[N-2] = (d[N-2]-b[N-2]*e[N-3])/c[N-2];

	b[N-1]  = b[N-1]-a[N-1]*d[N-3];
	c[N-1] = c[N-1]-a[N-1]*e[N-3]-b[N-1]*d[N-2];
}

// for NRHS = 1
void solve_penta(int N, double *alpha, double *beta, double *epsilon, 
				 double *gamma, double *delta, double *b){

	int i;

	// compute inverse(L) * b
	b[0] = b[0]/alpha[0];
	b[1] = (b[1]-beta[1]*b[0])/alpha[1];
	for (i=2; i<N; ++i){
		b[i] = (b[i]-epsilon[i]*b[i-2]-beta[i]*b[i-1])/alpha[i];
	}

	// compute inverse(U) * inverse(L) * b
	b[N-1] = b[N-1];
	b[N-2] = b[N-2]-gamma[N-2]*b[N-1];
	for (i=3; i<N+1; ++i){
		b[N-i] = b[N-i]-gamma[N-i]*b[N-i+1]-delta[N-i]*b[N-i+2];
	}
}

// NRHS >= 1
// b must be transpose beforehand
void solve_pentaX(int N, int NRHS, double *alpha, double *beta, double *epsilon, 
				  double *gamma, double *delta, double *b){

	int i, j;
	double *tempB;

	// solve column by column
	for (i=0; i<NRHS; ++i){

		// prepare RHS
		tempB = &b[i*N];
		// solve
		solve_penta(N, alpha, beta, epsilon, gamma, delta, tempB);
	}
}
