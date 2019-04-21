
#include <stdio.h>

void print_matrix(char *name, int N, int M, double *A){

    int i, j;

    printf("%s = \n", name);

    for (i=0; i<N; ++i){
            for (j=0; j<M; j++){

            	printf("%10f", A[i*M+j]);
            }
            printf( "\n" );
    }
    printf( "\n" );
}
