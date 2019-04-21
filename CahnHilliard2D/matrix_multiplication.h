
void matrix_multiplication(int M, int N, int K, double *A, double *B, double *C){

	// in : A is M*N
	// in : B is N*K
	// out: C is M*K

	int i, j, k;

	for (i=0; i<M; ++i){

		for (j=0; j<K; ++j){

			C[i*K+j] = 0;
			
			for (k=0; k<N; ++k){

				C[i*K+j] += A[i*N+k]*B[k*K+j];
			}
		}
	}
}
