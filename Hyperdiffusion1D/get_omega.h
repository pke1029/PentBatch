
void get_omega(double a, double b, double c, double d, double e, int N, double *E_inv, double *omega){

	// compute relavent entries of transpose(g) * (E^-1) * f
	double *gE;
	gE = (double*) calloc(2*4, sizeof(double));
	
	gE[0] = e*E_inv[0*N+0]     + a*E_inv[(N-2)*N+0]     + b*E_inv[(N-1)*N+0];
	gE[1] = e*E_inv[0*N+1]     + a*E_inv[(N-2)*N+1]     + b*E_inv[(N-1)*N+1];
	gE[2] = e*E_inv[0*N+(N-2)] + a*E_inv[(N-2)*N+(N-2)] + b*E_inv[(N-1)*N+(N-2)];
	gE[3] = e*E_inv[0*N+(N-1)] + a*E_inv[(N-2)*N+(N-1)] + b*E_inv[(N-1)*N+(N-1)];
	gE[4] = d*E_inv[0*N+0]     + e*E_inv[1*N+0]         + a*E_inv[(N-1)*N+0];
	gE[5] = d*E_inv[0*N+1]     + e*E_inv[1*N+1]         + a*E_inv[(N-1)*N+1];
	gE[6] = d*E_inv[0*N+(N-2)] + e*E_inv[1*N+(N-2)]     + a*E_inv[(N-1)*N+(N-2)];
	gE[7] = d*E_inv[0*N+(N-1)] + e*E_inv[1*N+(N-1)]     + a*E_inv[(N-1)*N+(N-1)];

	omega[0] = c - (a*gE[0] + e*gE[2] + d*gE[3]);
	omega[1] = d - (b*gE[0] + a*gE[1] + e*gE[3]);
	omega[2] = b - (a*gE[4] + e*gE[6] + d*gE[7]);
	omega[3] = c - (b*gE[4] + a*gE[5] + e*gE[7]);

	// invert omega
	double temp;
	double det = omega[0]*omega[3] - omega[1]*omega[2];
	temp = omega[0];
	omega[0] = omega[3] / det;
	omega[3] = temp / det;
	omega[1] = -omega[1] / det;
	omega[2] = -omega[2] / det;

	free(gE);
}
