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
ref_max_id = argmax(mat(ref_sig), rows=FALSE)
ref_min_id = argmin(mat(ref_sig), rows=FALSE)
sample_id[sample_num+1] = ref_max_id
sample_id[sample_num+2] = ref_min_id

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

### get scatterplot plot
tar_all = log2(tar_sig+small_num)
ref_all = log2(ref_sig+small_num)
tar_all_s = tar_all[sample_id]
ref_all_s = ref_all[sample_id]
tar_cpk_s = log2(tar_cpk_s+small_num)
ref_cpk_s = log2(ref_cpk_s+small_num)
tar_cbg_s = log2(tar_cbg_s+small_num)
ref_cbg_s = log2(ref_cbg_s+small_num)

all_upperlim = max(max(tar_all_s), max(ref_all_s))
all_lowerlim = min(min(tar_all_s), min(ref_all_s))

print(head(as.numeric(ref_all_s)))
print(head(as.numeric(tar_all_s)))
print(dim(cbind(ref_all_s, tar_all_s)))
print(all_lowerlim)
print(all_upperlim)
print(summary(cbind(ref_all_s, tar_all_s)))
png(paste(output, '.scatterplot.png', sep=''), width=800, height=400)
par(mfrow=c(1,2))
###### plot cluster
plot_MA_3parts(ref_all_s, tar_all_s, ref_cpk_s, tar_cpk_s, ref_cbg_s, tar_cbg_s, all_lowerlim, all_upperlim)
###### plot density plot
heatscatter(ref_all_s, tar_all_s, pch=16, ylim=c(all_lowerlim, all_upperlim), xlim=c(all_lowerlim, all_upperlim), cex=0.6)
abline(0, 1, col='black', lwd=2)
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
cpk_cbg_id_s = (cpk_id_s | cbg_id_s)

print(head(cbind(ref_all_s[cpk_cbg_id_s], tar_all_s[cpk_cbg_id_s])))

print(dim(cbind(ref_all_s[cpk_cbg_id_s], tar_all_s[cpk_cbg_id_s])))
print('get density weight (raw)')
density_weight_s = get_density_weight_2d(ref_all_s[cpk_cbg_id_s], tar_all_s[cpk_cbg_id_s], grid)



print('fit weighted loess model')
weighted_loess_model = loess(tar_all_s[cpk_cbg_id_s]~ref_all_s[cpk_cbg_id_s], weights=density_weight_s, span=1, degree=2)
print('fit loess model')
tar_all_pred = predict(weighted_loess_model, newdata=ref_all)
tar_all_pred_s = tar_all_pred[sample_id]
tar_all_norm1_s = tar_all_s + (ref_all_s-tar_all_pred_s)
tar_cpk_norm1_s = tar_cpk_s + (ref_all_s-tar_all_pred_s)[cpk_id_s]
tar_cbg_norm1_s = tar_cbg_s + (ref_all_s-tar_all_pred_s)[cbg_id_s]

print('get density weight (norm1)')
density_weight_norm1_s = get_density_weight_2d(ref_all_s[cpk_cbg_id_s], tar_all_norm1_s[cpk_cbg_id_s], grid)
#density_weight_norm1 = density_weight


png(paste(output, '.scatterplot.norm.png', sep=''), width=800, height=400)
par(mfrow=c(1,2))

###### plot cluster
#plot_MA_3parts(A_all, M_all, A_cpk, M_cpk, A_cbg, M_cbg, A_lim_lower, A_lim_upper)
heatscatter(ref_all_s, tar_all_s, pch=16, ylim=c(all_lowerlim, all_upperlim), xlim=c(all_lowerlim, all_upperlim), cex=0.6)
abline(0,1, col='black', lwd=2)
cpk_id_s_weight = cpk_id_s[cpk_cbg_id_s]
cbg_id_s_weight = cbg_id_s[cpk_cbg_id_s]
points(weighted_loess_model$x, weighted_loess_model$fitted,col="dodgerblue", pch=16, cex=0.6)
lines(c(weighted.mean(ref_cbg_s, density_weight_s[cbg_id_s_weight]), weighted.mean(ref_cpk_s, density_weight_s[cpk_id_s_weight])), c(weighted.mean(tar_cbg_s, density_weight_s[cbg_id_s_weight]), weighted.mean(tar_cpk_s, density_weight_s[cpk_id_s_weight])), col='green', lty=2, lwd=3)
points(weighted.mean(ref_cpk_s, density_weight_s[cpk_id_s_weight]), weighted.mean(tar_cpk_s, density_weight_s[cpk_id_s_weight]), pch=16, col='black', cex=1)
points(weighted.mean(ref_cbg_s, density_weight_s[cbg_id_s_weight]), weighted.mean(tar_cbg_s, density_weight_s[cbg_id_s_weight]), pch=16, col='black', cex=1)

###### plot loess adjusted MA-plot
heatscatter(ref_all_s, (tar_all_s+(ref_all_s-tar_all_pred_s)), pch=16, ylim=c(all_lowerlim, all_upperlim), xlim=c(all_lowerlim, all_upperlim), cex=0.6)
abline(0,1, col='black', lwd=2)
lines(c(weighted.mean(ref_cbg_s, density_weight_norm1_s[cbg_id_s_weight]), weighted.mean(ref_cpk_s, density_weight_norm1_s[cpk_id_s_weight])), c(weighted.mean(tar_cbg_norm1_s, density_weight_norm1_s[cbg_id_s_weight]), weighted.mean(tar_cpk_norm1_s, density_weight_norm1_s[cpk_id_s_weight])), col='green', lty=2, lwd=3)
points(weighted.mean(ref_cpk_s, density_weight_norm1_s[cpk_id_s_weight]), weighted.mean(tar_cpk_norm1_s, density_weight_norm1_s[cpk_id_s_weight]), pch=16, col='black', cex=1)
points(weighted.mean(ref_cbg_s, density_weight_norm1_s[cbg_id_s_weight]), weighted.mean(tar_cbg_norm1_s, density_weight_norm1_s[cbg_id_s_weight]), pch=16, col='black', cex=1)

dev.off()

###### get tar norm1 signal
#tar_sig_norm = (tar_sig) * (2^(ref_all-tar_all_pred))
tar_sig_norm = 2^(tar_all + ref_all - tar_all_pred) - small_num
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


