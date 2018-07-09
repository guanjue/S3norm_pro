###### NB random variable
nb_size = 0.9
nb_prob = 0.1
nb_n = 1e+6
nb_v = rnbinom(nb_n, nb_size, nb_prob)
nb_mean = mean(nb_v)
nb_var = var(nb_v)

#hist(nb_v, breaks=25)

zero_prob = 0.3
###### binary binomial variable
print(zero_prob)
binom_n = as.integer(nb_n * zero_prob / (1-zero_prob))
binom_v = rep(0.0, binom_n)
###### merge NB & binary
nb_binom_v = c(nb_v, binom_v)
###### plot histogram
#hist(nb_binom_v, breaks=25)
###### non-zero parameters
mean_all = mean(nb_binom_v)
var_all = var(nb_binom_v)
mean_non0 = mean(nb_binom_v[nb_binom_v>0])
var_non0 = var(nb_binom_v[nb_binom_v>0])
mean_x2_non0 = mean((nb_binom_v[nb_binom_v>0])^2)
###### non-zero Prob & Size
nb_prob_non0 = mean_non0 / var_non0
nb_size_non0 = mean_non0 * nb_prob_non0 / (1-nb_prob_non0)

###### identify p0
best_p0 = 0
best_prob_dif = 1

k=0
for (i in seq(0,0.99,0.005)){
	k = k+1
	ProbT = mean_non0 / (mean_x2_non0 - mean_non0^2 * (1-i))
	SizeT = mean_non0 * (1-i) * ProbT / (1-ProbT)
	nb_v_T = rnbinom(1e+4, SizeT, ProbT)
	p0_new = sum(nb_v_T==0) / length(nb_v_T)
	p0_dif = abs(i-p0_new)
	if ((k%%100)==0){
		print(paste('iteration:', toString(k)))
	}
	if (abs(i-p0_new) < best_prob_dif){
		print(paste('iteration:', toString(k)))
		print('change best_p0')
		best_prob_dif = abs(i-p0_new)
		best_p0 = i
	}
}

print('true p0:')
print(sum(nb_v==0)/length(nb_v))

print('estimated p0: ')
print(best_p0)



