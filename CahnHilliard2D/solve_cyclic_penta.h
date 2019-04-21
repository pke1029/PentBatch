
// solve cyclic pentadiagonal matrix given omega and LU factorisation of E
// temporary storage: d_hat(N-2), small_d(2)
// out: dn

void solve_cyclic_penta(int N, double *alpha, double *beta, double *epsilon, double *gamma, double *delta, 
						double *f, double *gE, double *omega, double *dn, double *small_d){
	
	int i;

	// compute transpose(g) * inverse(E) * d_hat
	matrix_multiplication(2, N-2, 1, gE, dn, small_d);

	// compute small_d - transpose(g) * inverse(E) * d_hat
	small_d[0] = dn[N-2] - small_d[0];
	small_d[1] = dn[N-1] - small_d[1];

	// compute omega * small_d
	dn[N-2] = omega[0]*small_d[0] + omega[1]*small_d[1];
	dn[N-1] = omega[2]*small_d[0] + omega[3]*small_d[1];

	// compute d_hat - f * small_x
	dn[0]   = dn[0]   - (f[0]*dn[N-2] + f[1]*dn[N-1]);
	dn[1]   = dn[1]   - f[3]*dn[N-1];
	dn[N-4] = dn[N-4] - f[(N-2)*2-4]*dn[N-2];
	dn[N-3] = dn[N-3] - (f[(N-2)*2-2]*dn[N-2] + f[(N-2)*2-1]*dn[N-1]);

	// compute inverse(E) * ( d_hat - * f * small_x )
	solve_penta(N-2, alpha, beta, epsilon, gamma, delta, dn);
}