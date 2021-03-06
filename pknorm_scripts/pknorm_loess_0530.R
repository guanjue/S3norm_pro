library(LSD)
library(MASS)
library(ramify)
library(mclust)

### get parameters
args = commandArgs(trailingOnly=TRUE)

ref = args[1]
tar = args[2]
output = args[3]

sample_num = as.numeric(args[4])
p_method = args[5]
fdr_thresh = as.numeric(args[6])
rank_top = as.numeric(args[7])

script_folder = args[8]


source(paste(script_folder, 'nbp_0326_0525.R', sep=''))

set.seed(2018)
### read signal
ref_sig = scan(ref)
tar_sig = scan(tar)

set.seed(2018)
sample_id = sample(length(ref_sig), sample_num)
A_max_id = argmax(mat(ref_sig+tar_sig), rows=FALSE)
A_min_id = argmin(mat(ref_sig+tar_sig), rows=FALSE)
sample_id[sample_num+1] = A_max_id
sample_id[sample_num+2] = A_min_id

### get nbp & pk
if (p_method == 'nb'){
	ref_p = get_2r_nbp(ref_sig, paste(ref, '.nbp.txt', sep=''))
	tar_p = get_2r_nbp(tar_sig, paste(tar, '.nbp.txt', sep=''))
	ref_p_pk_binary = p.adjust(ref_p, method='fdr') < fdr_thresh
	tar_p_pk_binary = p.adjust(tar_p, method='fdr') < fdr_thresh
}

ref_pk_num = sum(ref_p_pk_binary==1.0)
tar_pk_num = sum(tar_p_pk_binary==1.0)
print('ref_pk_num')
print(ref_pk_num)
print('tar_pk_num')
print(tar_pk_num)
### if smaller than 10000 pk, replace by top 10000
if (ref_pk_num<rank_top){
	print('ref: use top rank for peak')
	ref_p_pk_binary = ref_p<=sort(ref_p,decreasing=FALSE)[rank_top]
	ref_pk_num = sum(ref_p_pk_binary==1.0)
}
if (tar_pk_num<rank_top){
	print('tar: use top rank for peak')
	tar_p_pk_binary = ref_p<=sort(tar_p,decreasing=FALSE)[rank_top]
	tar_pk_num = sum(tar_p_pk_binary==1.0)
}

print('ref_pk_num')
print(ref_pk_num)
print('tar_pk_num')
print(tar_pk_num)
print('Adjusted random index')
ari = adjustedRandIndex(ref_p_pk_binary, tar_p_pk_binary)
print(ari)

### get common pk & bg ID
cpk_id = ref_p_pk_binary & tar_p_pk_binary
cbg_id = (ref_p_pk_binary + tar_p_pk_binary) == 0 
cpk_id_s = cpk_id[sample_id]
cbg_id_s = cbg_id[sample_id]

### get common pk & bg signal
ref_cpk_s = ref_sig[sample_id][cpk_id_s]
tar_cpk_s = tar_sig[sample_id][cpk_id_s]
ref_cbg_s = ref_sig[sample_id][cbg_id_s]
tar_cbg_s = tar_sig[sample_id][cbg_id_s]

ref_cpk_mean = mean(ref_cpk_s)
tar_cpk_mean = mean(tar_cpk_s)
ref_cbg_mean = mean(ref_cbg_s)
tar_cbg_mean = mean(tar_cbg_s)

### get added small number
small_num = (ref_cpk_mean*tar_cbg_mean - ref_cbg_mean*tar_cpk_mean) / ((ref_cbg_mean-ref_cpk_mean)-(tar_cbg_mean-tar_cpk_mean))
if (small_num >1){
	small_num = 1.0
} else if (small_num <0.01){
	small_num = 0.01
}
print(paste('added small number: ', toString(small_num), sep=''))

### get MA plot
M_all = (log2(tar_sig+small_num) - log2(ref_sig+small_num))
A_all = (0.5 * (log2(tar_sig+small_num) + log2(ref_sig+small_num)))
M_cpk_s = (log2(tar_cpk_s+small_num) - log2(ref_cpk_s+small_num))
A_cpk_s = (0.5 * (log2(tar_cpk_s+small_num) + log2(ref_cpk_s+small_num)))
M_cbg_s = (log2(tar_cbg_s+small_num) - log2(ref_cbg_s+small_num))
A_cbg_s = (0.5 * (log2(tar_cbg_s+small_num) + log2(ref_cbg_s+small_num)))

M_all_s = M_all[sample_id]
A_all_s = A_all[sample_id]


M_ylim = max(abs(min(M_all_s)), abs(max(M_all_s)))
A_lim_lower = min(A_all_s)
A_lim_upper = max(A_all_s)
png(paste(output, '.MAplot.png', sep=''), width=800, height=400)
par(mfrow=c(1,2))
###### plot cluster
plot_MA_3parts(A_all_s, M_all_s, A_cpk_s, M_cpk_s, A_cbg_s, M_cbg_s, M_ylim, A_lim_lower, A_lim_upper)
###### plot density plot
heatscatter(A_all_s, M_all_s, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
dev.off()


todiscrete = function(t, tmin, tmax, bins) {
	erg = round((t - tmin)/(tmax - tmin) * bins + 0.5)
	erg = pmin(pmax(erg, 1), bins)
	return(erg)
}


grid=100
cpk_cbg_id_s = (cpk_id_s | cbg_id_s)


print('fit weighted loess model')
loess_model = loess(M_all_s[cpk_cbg_id_s]~A_all_s[cpk_cbg_id_s], span=1, degree=2)


print('fit loess model')
M_all_pred = predict(loess_model, newdata=A_all)
M_all_pred_s = M_all_pred[sample_id]
M_all_norm1_s = M_all_s - M_all_pred_s
M_cpk_norm1_s = M_cpk_s - M_all_pred_s[cpk_id_s]
M_cbg_norm1_s = M_cbg_s - M_all_pred_s[cbg_id_s]

print('get density weight (norm1)')
density_weight_norm1_s = get_density_weight_2d(A_all_s[cpk_cbg_id_s], M_all_norm1_s[cpk_cbg_id_s], grid)
#density_weight_norm1 = density_weight


png(paste(output, '.MAplot.norm.png', sep=''), width=800, height=400)
par(mfrow=c(1,2))

###### plot cluster
#plot_MA_3parts(A_all, M_all, A_cpk, M_cpk, A_cbg, M_cbg, A_lim_lower, A_lim_upper)
heatscatter(A_all_s, M_all_s, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
cpk_id_s_weight = cpk_id_s[cpk_cbg_id_s]
cbg_id_s_weight = cbg_id_s[cpk_cbg_id_s]
points(loess_model$x, loess_model$fitted,col="dodgerblue", pch=16, cex=0.6)
lines(c(mean(A_cbg_s), mean(A_cpk_s)), c(mean(M_cbg_s), mean(M_cpk_s)), col='green', lty=2, lwd=3)
points(mean(A_cpk_s), mean(M_cpk_s), pch=16, col='black', cex=1)
points(mean(A_cbg_s), mean(M_cbg_s), pch=16, col='black', cex=1)

###### plot loess adjusted MA-plot
heatscatter(A_all_s, (M_all_s-M_all_pred_s), pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), cex=0.6)
abline(h=0, col='black', lwd=2)
lines(c(mean(A_cbg_s), mean(A_cpk_s)), c(mean(M_cbg_norm1_s), mean(M_cpk_norm1_s)), col='green', lty=2, lwd=3)
points(mean(A_cpk_s), mean(M_cpk_norm1_s), pch=16, col='black', cex=1)
points(mean(A_cbg_s), mean(M_cbg_norm1_s), pch=16, col='black', cex=1)

dev.off()

###### get tar norm1 signal
tar_sig_norm = (tar_sig) / (2^M_all_pred)

write.table(tar_sig_norm, output, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### get nbp & pk
if (p_method == 'nb'){
	tar_norm_p = get_2r_nbp(tar_sig_norm, paste(tar, '.norm.nbp.txt', sep=''))
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


