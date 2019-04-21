
void get_penta_mat(double a, double b, double c, double d, double e, int N, double *A){

	int i, j;

	A[0*N+0] = c;
	A[0*N+1] = d;
	A[0*N+2] = e;

	A[1*N+0] = b;
	A[1*N+1] = c;
	A[1*N+2] = d;
	A[1*N+3] = e;

	for (i=2; i<N-2; ++i){
		
		A[i*N+(i-2)] = a;
		A[i*N+(i-1)] = b;
		A[i*N+(i-0)] = c;
		A[i*N+(i+1)] = d;
		A[i*N+(i+2)] = e;
	}

	A[(N-2)*N+(N-4)] = a;
	A[(N-2)*N+(N-3)] = b;
	A[(N-2)*N+(N-2)] = c;
	A[(N-2)*N+(N-1)] = d;

	A[(N-1)*N+(N-3)] = a;
	A[(N-1)*N+(N-2)] = b;
	A[(N-1)*N+(N-1)] = c;
}
