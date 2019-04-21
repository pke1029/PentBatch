
// compute omega for a cyclic pentadiagonal matrix
// alpha, beta, epsilon, gamma and delta are the LU fators of the dim N-2 (non-cyclic) pentadiagonal matrix

void get_omega(int N2, double b, double c, double d,
			   double *alpha, double *beta, double *epsilon, double *gamma, 
			   double *delta, double *gE, double *f, double *omega){

	double temp;

	// compute transpose(g) * inverse(E) * f
	matrix_multiplication(2, N2, 2, gE, f, omega);

	// compute W - transpose(g) * inverse(E) * f
	omega[0] = c - omega[0];
	omega[1] = d - omega[1];
	omega[2] = b - omega[2];
	omega[3] = c - omega[3];

	// invert omega
	double det = omega[0]*omega[3] - omega[1]*omega[2];
	temp = omega[0];
	omega[0] = omega[3] / det;
	omega[3] = temp / det;
	omega[1] = -omega[1] / det;
	omega[2] = -omega[2] / det;
}
