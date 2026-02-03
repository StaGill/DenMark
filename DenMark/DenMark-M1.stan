///////////////////////////////////////////////////////////////////
// Model 1: two independent latent fields, drop a21 field dependence
// log(lambda1) = beta0 + w1(g)
// log(lambda2) = beta1 + log(Y) +  a22*w2(g)
///////////////////////////////////////////////////////////////////


functions {
	// compute the eigenvalues: 
	vector compute_lambda(array[] real L, array[] int m, int d) {
		vector[d] lambda;
		for(i in 1:d){
			lambda[i] = ((m[i]*pi())/(2*L[i]))^2;
			}
		return lambda;
	}
	// matern 32 
	real spd_isotropic_matern32(real sigma, real ell, vector x, int d) {
	  real spd = square(sigma) * 6.0 * pi() * pow(3.0, 1.5) * pow(ell, d) * inv(pow(3 + square(ell) * dot_self(x), 2.5));
	  return spd;
	  }

	vector compute_PHI(array[] real L, array[] int m, matrix coords) {
		int c = cols(coords);
		int r = rows(coords);
		
		matrix[r,c] phi;
		vector[r] phi1;
		for (i in 1:c){
			phi[,i] = 1/sqrt(L[i])*sin(m[i]*pi()*(coords[,i]+L[i])/(2*L[i]));
		}
		phi1 = phi[,1];
		for (i in 2:c){
			phi1 = phi1 .* phi[,i];
		}
		return phi1;
	}

}

data {
	int<lower=1> d;
	array[d] real L;
	int<lower=1> mstar;
	int<lower=1> n;
	matrix[n,d] coords;
	real log_grid_area;
	array[n] int y1;
	array[n] int y2;
	real<lower=0.00001>mrange;
	
	array[mstar,d] int indices;
	int is_centerted_PHI;
}

transformed data {
	matrix[n,mstar] PHI;
	for (j in 1:mstar){ 
	  PHI[,j] = compute_PHI(L, indices[j,], coords); 
	}
	if(is_centerted_PHI==1){
	  matrix[n,mstar] cPHI;
	  for(i in 1:n){
	    cPHI[i,] = PHI[i,] - mean(PHI[i,]);
	}
	PHI = cPHI;
}
}

parameters {
  vector<lower=0.0001>[2] ell;
	real beta0;
	real beta1;
	real<lower=0.0001> a11;
	real<lower=0.0001> a22;
	vector[mstar] betab1;
	vector[mstar] betab2;
}

transformed parameters{
	vector[mstar] diagSPD1;
	vector[mstar] diagSPD2;
	for(i in 1:mstar){ 
	  diagSPD1[i] =  sqrt(spd_isotropic_matern32(1, ell[1], sqrt(compute_lambda(L, indices[i,], d)), d)); // (8)
	  diagSPD2[i] =  sqrt(spd_isotropic_matern32(1, ell[2], sqrt(compute_lambda(L, indices[i,], d)), d));
	}
	vector[n] w1 = PHI * (diagSPD1 .* betab1); // f(x)~ sum (f(x)*SPD), (8)
	vector[n] w2 = PHI * (diagSPD2 .* betab2);
}

model{
	beta0 ~ normal(0, 5); // intercept 1  
    beta1 ~ normal(0, 5); // intercept 2 
    a11 ~ normal(0, 0.5); // SD of the 1st GP
    a22 ~ normal(0, 0.5); // SD of the 2nd GP
	ell ~ inv_gamma(2, mrange); // length of scale for two GPs
	
	betab1 ~ std_normal();
	betab2 ~ std_normal();
	
	vector[n] z1 = a11*w1;
	vector[n] z2 = a22*w2;
	
	y1 ~ poisson_log( log_grid_area + beta0 + z1 );
	y2 ~ poisson_log( beta1 + log(to_vector(y1) + 0.05 ) + z2  );
}

generated quantities{
  vector[n] log_lik;
  {
    vector[n] z1 = a11*w1;
    vector[n] z2 = a22*w2;
    for (i in 1:n) {
      real lambda1 = log_grid_area + beta0 + z1[i];
      real lambda2 = beta1 + log(y1[i] + 0.05) + z2[i];
     // Joint log-likelihood for WAIC
      log_lik[i] = poisson_log_lpmf(y1[i] | lambda1) +
                   poisson_log_lpmf(y2[i] | lambda2);
    }
  }
}





