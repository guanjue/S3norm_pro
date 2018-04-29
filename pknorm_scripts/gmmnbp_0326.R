library(mixtools)
set.seed(2018)

### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
output_name = args[2]

### remove 0 bins
sig = scan(signal_track_file)
thesh = 0
sig_non0 = sig[sig>thesh]
### log2 transform
log_sig_non0 = log2(sig_non0)

### fit gmm 2c
cluster_num = 2
### random select 20% of data to fit gmm
sample_num = length(log_sig_non0)/5
sample_id = sample(length(log_sig_non0), sample_num)
log_sig_non0_sample = log_sig_non0[sample_id,]
### fit gmm
mixmdl = normalmixEM(log_sig_non0, k=cluster_num)
gmm_posterior = mixmdl$posterior
gmm_cluster = apply(gmm_posterior, 1, function(x) which(x==max(x)))
mean_c = mixmdl$mu
pdf(paste(output_name, '.hist_gmm.pdf', sep=''), width=6, height=6)
plot(mixmdl,which=2)
dev.off()

if (mean_c[1] <= mean_c[2]){
	sig_mean = mean_c[1]
	bg_bins = log_sig_non0_sample[gmm_cluster==1,]
	sig_var = var(bg_bins)
} else{
	sig_mean = mean_c[2]
	bg_bins = log_sig_non0_sample[gmm_cluster==2,]
	sig_var = var(bg_bins)
}

sig_mean_raw = mean(sig_non0)
sig_var_raw = var(sig_non0)
print(paste('check signal track overdispersion in background regions (raw), var/mean=', toString(round(sig_var_raw/sig_mean_raw, digits=3)) ))
print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_var/sig_mean, digits=3)) ))

sig_prob = sig_mean / sig_var
if (sig_prob<0.1){
	print('use lower bound for nbp prob')
	sig_prob = 0.1
}

if (sig_prob>=0.9){
	print('use upper bound for nbp prob')
	sig_prob = 0.9
}

sig_size = sig_mean * sig_prob / (1-sig_prob)

nb_pval = apply(as.matrix(sig), MARGIN=1, function(x) pnbinom(x[1], sig_size, sig_prob, lower.tail=FALSE) )

### write output
write.table(nb_pval, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

