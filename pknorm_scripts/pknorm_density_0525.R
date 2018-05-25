library(LSD)
library(MASS)
library(ramify)
library(mclust)

### get parameters
args = commandArgs(trailingOnly=TRUE)

ref = args[1]
tar = args[2]
output = args[3]

sample_num = args[4]
p_method = args[5]
fdr_thresh = args[6]

script_folder = args[7]

source(paste(script_folder, 'nbp_0326_0525.R', sep=''))

set.seed(2018)
### read signal
ref_sig = scan(ref)
tar_sig = scan(tar)

### get nbp & pk
if (p_method == 'nb'){
	ref_p = get_2r_nbp(ref_sig, paste(ref, '.nbp.txt', sep=''))
	tar_p = get_2r_nbp(tar_sig, paste(ref, '.nbp.txt', sep=''))
	ref_p_pk_binary = p.adjust(ref_p, method='fdr') < fdr_thresh
	tar_p_pk_binary = p.adjust(tar_p, method='fdr') < fdr_thresh
}

ref_pk_num = sum(ref_p_pk_binary==1.0)
print('ref_pk_num')
print(ref_pk_num)
tar_pk_num = sum(tar_p_pk_binary==1.0)
print('tar_pk_num')
print(tar_pk_num)
print('Adjusted random index')
ari = adjustedRandIndex(ref_p_pk_binary, tar_p_pk_binary)
print(ari)

### get common pk & bg ID
cpk_id = ref_p_pk_binary & tar_p_pk_binary
cbg_id = (ref_p_pk_binary + tar_p_pk_binary) == 0 

### get common pk & bg signal
ref_cpk = ref_sig[cpk_id]
tar_cpk = tar_sig[cpk_id]
ref_cbg = ref_sig[cbg_id]
tar_cbg = tar_sig[cbg_id]

ref_cpk_mean = mean(ref_cpk)
tar_cpk_mean = mean(tar_cpk)
ref_cbg_mean = mean(ref_cbg)
tar_cbg_mean = mean(tar_cbg)

### get added small number
small_num = (ref_cpk_mean*tar_cbg_mean - ref_cbg_mean*tar_cpk_mean) / ((ref_cbg_mean-ref_cpk_mean)-(tar_cbg_mean-tar_cpk_mean))
if (small_num >1){
	small_num = 1.0
} else if (small_num <0.01){
	small_num = 0.01
}
print(paste('added small number: ', toString(small_num), sep=''))

### get MA plot
M_all = log2(tar_sig+small_num) - log2(ref_sig+small_num)
A_all = 0.5 * (log2(tar_sig+small_num) + log2(ref_sig+small_num))
M_cpk = log2(tar_cpk+small_num) - log2(ref_cpk+small_num)
A_cpk = 0.5 * (log2(tar_cpk+small_num) + log2(ref_cpk+small_num))
M_cbg = log2(tar_cbg+small_num) - log2(ref_cbg+small_num)
A_cbg = 0.5 * (log2(tar_cbg+small_num) + log2(ref_cbg+small_num))

M_ylim = max(abs(min(M_all)), abs(max(M_all)))
A_lim_lower = min(A_all)
A_lim_upper = max(A_all)
png(paste(tar, '.MAplot.png', sep=''), width=800, height=400)
par(mfrow=c(1,2))
###### plot cluster
plot_MA_3parts(A_all, M_all, A_cpk, M_cpk, A_cbg, M_cbg, A_lim_lower, A_lim_upper)
###### plot density plot
heatscatter(A_all, M_all, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
dev.off()


todiscrete = function(t, tmin, tmax, bins) {
	erg = round((t - tmin)/(tmax - tmin) * bins + 0.5)
	erg = pmin(pmax(erg, 1), bins)
	return(erg)
}


getfrommat_2d = function(d, a) {
	d$z[a[1], a[2]]
}

getfrommat_1d = function(d, a) {
	d$y[a[1]]
}

get_density_weight_2d = function(A_all, M_all, grid){
	###### discrete original 2D space
	A_all_discrete = todiscrete(A_all, min(A_all), max(A_all), bins = grid)
	M_all_discrete = todiscrete(M_all, min(M_all), max(M_all), bins = grid)
	d = kde2d(A_all, M_all, n = grid)
	###### get density
	density = unlist(apply(cbind(A_all_discrete, M_all_discrete), 1, function(x) getfrommat_2d(d,x)))
	density_weight_raw = 1/(density)
	#density_weight_raw[cpk_id] = mean(density_weight_raw[cpk_id])
	#density_weight_raw = density_weight_raw - min(density_weight_raw)
	density_weight = density_weight_raw/sum(density_weight_raw)
	return(density_weight)
}

get_density_weight_1d = function(A_all, grid){
	###### discrete original 2D space
	A_all_discrete = todiscrete(A_all, min(A_all), max(A_all), bins = grid)
	d = density(A_all, n = grid)
	###### get density
	density = unlist(lapply(A_all_discrete, function(x) getfrommat_1d(d,x)))
	density_weight_raw = 1/(density)
	#density_weight_raw[cpk_id] = mean(density_weight_raw[cpk_id])
	#density_weight_raw = density_weight_raw - min(density_weight_raw)
	density_weight = density_weight_raw/sum(density_weight_raw)
	return(density_weight)
}

grid=100
cpk_cbg_id = cpk_id | cbg_id

print('get density weight (raw)')
density_weight = get_density_weight_2d(A_all, M_all, grid)


set.seed(2018)
sample_id = sample(length(cpk_cbg_id), sample_num)
A_max_id = argmax(mat(A_all[cpk_cbg_id]), rows=FALSE)
A_min_id = argmin(mat(A_all[cpk_cbg_id]), rows=FALSE)
sample_id[sample_num+1] = A_max_id
sample_id[sample_num+2] = A_min_id

print('fit weighted loess model')
weighted_loess_model = loess(M_all[cpk_cbg_id][sample_id]~A_all[cpk_cbg_id][sample_id], weights=density_weight[cpk_cbg_id][sample_id], span=1, degree=2)
print('fit loess model')
loess_model = loess(M_all[sample_id]~A_all[sample_id], span=1, degree=2, control=loess.control(surface="direct"))
M_all_pred = predict(weighted_loess_model, newdata=A_all)
M_all_norm1 = M_all - M_all_pred
M_cpk_norm1 = M_cpk - M_all_pred[cpk_id]
M_cbg_norm1 = M_cbg - M_all_pred[cbg_id]

print('get density weight (norm1)')
density_weight_norm1 = get_density_weight_2d(A_all, M_all_norm1, grid)
#density_weight_norm1 = density_weight


png(paste(tar, '.MAplot.norm.png', sep=''), width=800, height=800)
par(mfrow=c(2,2))

###### plot cluster
#plot_MA_3parts(A_all, M_all, A_cpk, M_cpk, A_cbg, M_cbg, A_lim_lower, A_lim_upper)
heatscatter(A_all, M_all, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
points(weighted_loess_model$x, weighted_loess_model$fitted,col="dodgerblue", pch=16, cex=0.6)
lines(c(weighted.mean(A_cbg, density_weight[cbg_id]), weighted.mean(A_cpk, density_weight[cpk_id])), c(weighted.mean(M_cbg, density_weight[cbg_id]), weighted.mean(M_cpk, density_weight[cpk_id])), col='green', lty=2, lwd=3)
points(weighted.mean(A_cpk, density_weight[cpk_id]), weighted.mean(M_cpk, density_weight[cpk_id]), pch=16, col='black', cex=1)
points(weighted.mean(A_cbg, density_weight[cbg_id]), weighted.mean(M_cbg, density_weight[cbg_id]), pch=16, col='black', cex=1)

#plot_MA_3parts(A_all, M_all, A_cpk, M_cpk, A_cbg, M_cbg, A_lim_lower, A_lim_upper)
heatscatter(A_all, M_all, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
points(as.vector(loess_model$x), as.vector(loess_model$fitted),col="dodgerblue", pch=16, cex=0.6)
lines(c(mean(A_cbg), mean(A_cpk)), c(mean(M_cbg), mean(M_cpk)), col='green', lty=2, lwd=3)
points(mean(A_cpk), mean(M_cpk), pch=16, col='black', cex=1)
points(mean(A_cbg), mean(M_cbg), pch=16, col='black', cex=1)

###### plot loess adjusted MA-plot
heatscatter(A_all, M_all-M_all_pred, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
lines(c(weighted.mean(A_cbg, density_weight_norm1[cbg_id]), weighted.mean(A_cpk, density_weight_norm1[cpk_id])), c(weighted.mean(M_cbg_norm1, density_weight_norm1[cbg_id]), weighted.mean(M_cpk_norm1, density_weight_norm1[cpk_id])), col='green', lty=2, lwd=3)
points(weighted.mean(A_cpk, density_weight_norm1[cpk_id]), weighted.mean(M_cpk_norm1, density_weight_norm1[cpk_id]), pch=16, col='black', cex=1)
points(weighted.mean(A_cbg, density_weight_norm1[cbg_id]), weighted.mean(M_cbg_norm1, density_weight_norm1[cbg_id]), pch=16, col='black', cex=1)

dev.off()

###### get tar norm1 signal
tar_sig_norm = (tar_sig) / (2^M_all_pred)

write.table(tar_sig_norm, output, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### get nbp & pk
if (p_method == 'nb'){
	tar_norm_p = get_2r_nbp(tar_sig_norm, paste(ref, '.nbp.txt', sep=''))
	tar_norm_p_pk_binary = p.adjust(tar_norm_p, method='fdr') < fdr_thresh
}

ref_pk_num = sum(ref_p_pk_binary==1.0)
print('ref_pk_num')
print(ref_pk_num)
tar_pk_num = sum(tar_p_pk_binary==1.0)
print('tar_pk_num')
print(tar_pk_num)
tar_norm_pk_num = sum(tar_norm_p_pk_binary==1.0)
print('tar_norm_pk_num')
print(tar_norm_pk_num)
print('Adjusted random index')
ari = adjustedRandIndex(ref_p_pk_binary, tar_norm_p_pk_binary)
print(ari)

#time Rscript /Users/universe/Documents/2018_BG/PKnorm/pknorm_scripts/pknorm_density_0525.R t.1.txt t.2.txt t.2.txt.pknorm.txt
#time Rscript /Users/universe/Documents/2018_BG/PKnorm/pknorm_scripts/pknorm_density_0525.R t.1.txt t.2.txt.pknorm.txt t.2.txt.pknorm.txt.pknorm.txt


