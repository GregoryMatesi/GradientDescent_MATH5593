param N >= 0;	# 20 genetic allele frequency observations
set ancestry;	# 5 global ancestries in the mixtures model

param a_tilde {1..N};	# 20x1 allele frequency from an observed heterogeneous ancestry
param A{1..N, ancestry};			# 20x5 matrix of allele frequencies

var pi {ancestry} >= 0;	# propotion value for each ancestry in the mixtures model

minimize HA:
	sum {n in 1..N}  (sum {k in ancestry} ((A[n,k] * pi[k]) - a_tilde[n])**2) ;
	# sum of least squares mixtures model
	
subject to sumtoone:
	sum {k in ancestry} pi[k] = 1;	# proportions sum to 1
