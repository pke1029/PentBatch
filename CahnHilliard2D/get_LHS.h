
void get_penta_mat(int N2, double a, double b, double c, double d, double e, 
				   double *a_vec, double *b_vec, double *c_vec, double *d_vec, double *e_vec){

	int i;

	c_vec[0] = c;
	d_vec[0] = d;
	e_vec[0] = e;

	b_vec[1] = b;
	c_vec[1] = c;
	d_vec[1] = d;
	e_vec[1] = e;

	for (i=2; i<N2-2; ++i){

		a_vec[i] = a;
		b_vec[i] = b;
		c_vec[i] = c;
		d_vec[i] = d;
		e_vec[i] = e;
	}

	a_vec[N2-2] = a;
	b_vec[N2-2] = b;
	c_vec[N2-2] = c;
	d_vec[N2-2] = d;

	a_vec[N2-1] = a;
	b_vec[N2-1] = b;
	c_vec[N2-1] = c;
}

void get_fg(int N2, double a, double b, double d, double e, double *f, double *g){

	int i;
	for (i=0; i<N2*2; ++i){
		f[i] = 0;
		g[i] = 0;
	}

	f[0]      = a;
	f[1]      = b;
	f[3]      = a;
	f[2*N2-4] = e;
	f[2*N2-2] = d;
	f[2*N2-1] = e;

	g[0] 	  = e;
	g[N2-2]   = a;
	g[N2-1]   = b;
	g[N2] 	  = d;
	g[N2+1]   = e;
	g[2*N2-1] = a;
}
