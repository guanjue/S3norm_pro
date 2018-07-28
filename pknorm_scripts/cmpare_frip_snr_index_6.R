library(mclust)
#library(Hmisc)

### get parameters
args = commandArgs(trailingOnly=TRUE)

sig1_raw = args[1]
sig2_raw = args[2]

sig2_TRnorm = args[3]

sig2_MAnorm = args[4]

sig2_QTnorm = args[5]

sig2_PKnorm = args[6]

sig1_Znorm = args[7]
sig2_Znorm = args[8]

output_name = args[9]
fdr_thresh = as.numeric(args[10])
method=args[11]

bed_file=args[12]

upplim = as.numeric(args[13])

#sig1_raw = 'ERY_ad.ctcfrep.fisher_p.txt'
#sig2_raw = 'T_CD8_SPL.ctcfrep.fisher_p.txt'

#sig2_TRnorm = 'T_CD8_SPL.trnorm.txt'

#sig2_MAnorm = 'T_CD8_SPL.manorm.txt'

#sig2_QTnorm = 'T_CD8_SPL.ctcfrep.fisher_p.txt.qtn.txt'

#sig2_PKnorm = 'T_CD8_SPL_fisher_p.pknorm.txt'

#sig1_Znorm = 'ERY_ad.znorm.txt'
#sig2_Znorm = 'T_CD8_SPL.znorm.txt'

#output_name = 'T_CD8_SPL.ctcfrep.fisher_p.txt.all5info.txt'
#fdr_thresh = as.numeric(0.05)
#method='neglog10p'

#bed_file='/storage/home/gzx103/scratch/vision/all_final_data/200_noblack.11_22_2017.bed'


################################################
nbp_2r = function(sig, p_lim_1r, output_name){
	### get mean & var
	thesh = -1
	sig_non0 = sig[sig>=thesh]
	sig_mean = mean(sig_non0)
	sig_var = var(sig_non0)
	### print overdispersion
	print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_var/sig_mean, digits=3)) ))
	### get negative binomial prob
	sig_prob = sig_mean / sig_var
	### set upper & lower lim for the NB-prob
	if (sig_prob<0.1){
		print('use lower bound for nbp prob')
		sig_prob = 0.1
	}
	if (sig_prob>=0.9){
		print('use upper bound for nbp prob')
		sig_prob = 0.9
	}
	### get size
	sig_size = sig_mean * sig_prob / (1-sig_prob)
	### get NB-p-value
	nb_pval = apply(as.matrix(sig_non0), MARGIN=1, function(x) pnbinom(x[1], sig_size, sig_prob, lower.tail=FALSE) )
	### 2nd round
	sig_non0_bg = sig_non0[nb_pval>=p_lim_1r]
	sig_non0_bg_mean = mean(sig_non0_bg)
	sig_non0_bg_var = var(sig_non0_bg)
	print(paste('2nd round check signal track overdispersion in background regions, var/mean=', toString(round(sig_var/sig_mean, digits=3)) ))
	sig_prob_2nd = sig_non0_bg_mean / sig_non0_bg_var
	if (sig_prob_2nd<0.1){
		print('use lower bound for nbp prob')
		sig_prob_2nd = 0.1
	}
	if (sig_prob_2nd>=0.9){
		print('use upper bound for nbp prob')
		sig_prob_2nd = 0.9
	}
	sig_size_2nd = sig_non0_bg_mean * sig_prob_2nd / (1-sig_prob_2nd)
	nb_pval_2nd = apply(as.matrix(sig), MARGIN=1, function(x) pnbinom(x[1], sig_size_2nd, sig_prob_2nd, lower.tail=FALSE) )

	### write nbp_2r
	#write.table(nb_pval_2nd, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
	### return 2nd round NB-p-value
	return(nb_pval_2nd)
}


zp = function(sig){
	### get mean & var
	thesh = -1
	sig_non0 = sig[sig>=thesh]
	sig_mean = mean(sig_non0)
	sig_sd = sd(sig_non0)
	z = (sig - sig_mean)/sig_sd
	zp = pnorm(-z)
	return(zp)
}


zp_2r = function(sig, p_thresh){
	### get mean & var
	thesh = -1
	sig_non0 = sig[sig>=thesh]
	sig_mean = mean(sig_non0)
	sig_sd = sd(sig_non0)
	z = (sig - sig_mean)/sig_sd
	zp = pnorm(-z)
	sig_2r = sig_non0[zp>=p_thresh]
	sig_mean = mean(sig_2r)
	sig_sd = sd(sig_2r)
	z = (sig - sig_mean)/sig_sd
	zp = pnorm(-z)	
	return(zp)
}


jaccard_index = function(pk_binary_1, pk_binary_2){
	overlap = (pk_binary_1 * pk_binary_2)==1
	pk_overlap_num = sum(overlap)
	union = (pk_binary_2 + pk_binary_2)!=0
	pk_union_num = sum(union)
	jaccard_index = pk_overlap_num / pk_union_num
	return(jaccard_index)
}


frip_common = function(sig1, pk_binary_1, sig2, pk_binary_2){
	overlap = (pk_binary_1 * pk_binary_2)==1
	union = (pk_binary_2 + pk_binary_2)!=0
	FRiP_JI_all = (sum(sig1[overlap])+sum(sig2[overlap])) / (sum(sig1)+sum(sig2))
	FRiP_JI_x = sum(sig1[overlap]) / sum(sig1)
	FRiP_JI_y = sum(sig2[overlap]) / sum(sig2)
	return(c(FRiP_JI_all, FRiP_JI_x, FRiP_JI_y))
}


plot_scatterplot_3parts = function(ref_all_s, tar_all_s, ref_cpk_s, tar_cpk_s, ref_cbg_s, tar_cbg_s, cpk_mean, cbg_mean, all_mean, all_lowerlim, all_upperlim){
	plot(ref_all_s, tar_all_s, pch=16, ylim=c(all_lowerlim, all_upperlim), xlim=c(all_lowerlim, all_upperlim), col='dodgerblue', cex=0.6)
	points(ref_cpk_s, tar_cpk_s, pch=16, col='darkorange1', cex=0.6)
	points(ref_cbg_s, tar_cbg_s, pch=16, col='gray56', cex=0.6)
	points(cpk_mean[1], cpk_mean[2], pch=16, col='black', cex=1)
	points(cbg_mean[1], cbg_mean[2], pch=16, col='black', cex=1)
	points(all_mean[1], all_mean[2], pch=16, col='red', cex=1.2)
	lines(c(cbg_mean[1], cpk_mean[1]), c(cbg_mean[2], cpk_mean[2]), col='green', lty=2, lwd=3)
	abline(0,1, col='black', lwd=3)
}

#sig1_raw = 'ER4.fisher_p.txt'
#sig2_raw = 'B_SPL.fisher_p.txt'
#sig2_TRnorm = 'B_SPL_TM_MAnorm_totalmean_MAnorm/B_SPL_TM_MAnorm.totalsig_norm.txt'
#sig2_MAnorm = 'B_SPL_TM_MAnorm_totalmean_MAnorm/B_SPL_TM_MAnorm.MAnorm.txt'
#sig2_QTnorm = 'B_SPL.fisher_p.txt.qtn.txt'
#sig2_PKnorm = 'B_SPL_txt.pknorm.txt'
#output_name = 'test'
#fdr_thresh = 0.05

### read data
sig1 = scan(sig1_raw)

sig2_r = scan(sig2_raw)
sig2_r[sig2_r>upplim] = upplim
sig2_tr = scan(sig2_TRnorm)
sig2_tr[sig2_tr>upplim] = upplim
sig2_ma = scan(sig2_MAnorm)
sig2_MAnorm[sig2_MAnorm>upplim] = upplim
sig2_qt = scan(sig2_QTnorm)
sig2_QTnorm[sig2_QTnorm>upplim] = upplim
sig2_pk = scan(sig2_PKnorm)
sig2_PKnorm[sig2_PKnorm>upplim] = upplim
sig1_z = scan(sig1_Znorm)
#sig1_Znorm[sig1_Znorm>upplim] = upplim
sig2_z = scan(sig2_Znorm)
#sig2_Znorm[sig2_Znorm>upplim] = upplim

sig_matrix = as.matrix(cbind(sig1, sig2_r, sig2_tr, sig2_ma, sig2_qt, sig2_pk, sig1_z, sig2_z))
sig_matrix_colnames = c('sig1', 'sig2_r', 'sig2_tr', 'sig2_ma', 'sig2_qt', 'sig2_pk', 'sig1_z', 'sig2_z')

### get nbp
if (method == 'nbp'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(nbp_2r(x, 0.001, paste(sig1_raw, '.nbp_2r.txt', sep='')), method='fdr'))
} else if (method == 'p'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(x, method='fdr'))	
} else if (method == 'neglog10p'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(10^(-x), method='fdr'))
} else if (method == 'zp'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(zp(x), method='fdr'))
} else if (method == 'zp_2r'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(zp_2r(x, 0.001), method='fdr'))
}

sig_matrix_p[,dim(sig_matrix_p)[2]-1] = p.adjust(pnorm(-sig1_z), method='fdr')
sig_matrix_p[,dim(sig_matrix_p)[2]] = p.adjust(pnorm(-sig2_z), method='fdr')



frip_all = c()
snr_all = c()
pk_num_all = c()
ari_all = c()
ji_all = c()

frip_common_all = c()

sig1_pk_id = sig_matrix_p[,1]<fdr_thresh
sig1 = sig_matrix[,1]
sig1_pk_all = c()
sig1_mean_pk = mean(sig1[sig1_pk_id])
sig1_mean_bg = mean(sig1[!sig1_pk_id])

set.seed(2018)
sample_id = sample(dim(sig_matrix)[1], 500000)

for ( i in c(1:(dim(sig_matrix)[2]-2))){
	print(sig_matrix_colnames[i])
	### get info
	nbp_tmp = sig_matrix_p[,i]
	sig_tmp = sig_matrix[,i]
	pk_id_tmp = nbp_tmp<fdr_thresh
	sig1_pk_all = cbind(sig1_pk_all, pk_id_tmp)
	bg_id_tmp = nbp_tmp>=fdr_thresh
	frip_tmp = sum(sig_tmp[pk_id_tmp]) / sum(sig_tmp)
	mean_pk = mean(sig_tmp[pk_id_tmp])
	mean_bg = mean(sig_tmp[bg_id_tmp])
	snr_tmp = mean_pk / mean_bg
	pk_num_tmp = sum(pk_id_tmp)
	ari_tmp = adjustedRandIndex(sig1_pk_id, pk_id_tmp)
	ji_tmp = jaccard_index(sig1_pk_id, pk_id_tmp)
	frip_common_tmp = frip_common(sig1, sig1_pk_id, sig_tmp, pk_id_tmp)
	### append(info)
	frip_all[i] = frip_tmp
	snr_all[i] = snr_tmp
	pk_num_all[i] = pk_num_tmp
	ari_all[i] = ari_tmp
	ji_all[i] = ji_tmp
	frip_common_all = cbind(frip_common_all, frip_common_tmp)
	### plotting
	cpk = (sig1_pk_id * pk_id_tmp) == 1
	cbg = (sig1_pk_id + pk_id_tmp) == 0
	ref_all_s = log10(sig1[sample_id]+0.1)
	tar_all_s = log10(sig_tmp[sample_id]+0.1)
	cpk_s = cpk[sample_id]
	cbg_s = cbg[sample_id]
	ref_cpk_s = ref_all_s[cpk_s]
	tar_cpk_s = tar_all_s[cpk_s]
	ref_cbg_s = ref_all_s[cbg_s]
	tar_cbg_s = tar_all_s[cbg_s]
	all_mean = c(log10(mean(sig1)+0.1), log10(mean(sig_tmp)+0.1))
	cpk_mean = c(log10(mean(sig1[cpk])+0.1), log10(mean(sig_tmp[cpk])+0.1))
	cbg_mean = c(log10(mean(sig1[cbg])+0.1), log10(mean(sig_tmp[cbg])+0.1))
	png(paste(output_name, sig_matrix_colnames[i], '.scatter.png', sep=''))
	plot_scatterplot_3parts(ref_all_s, tar_all_s, ref_cpk_s, tar_cpk_s, ref_cbg_s, tar_cbg_s, cpk_mean, cbg_mean, all_mean, -1.1, 2.5)
	dev.off()
}

### sig1z
nbp_tmp = sig_matrix_p[,dim(sig_matrix_p)[2]-1]
sig_tmp = sig_matrix[,dim(sig_matrix_p)[2]-1]
pk_id_tmp = nbp_tmp<fdr_thresh
sig1_pk_all = cbind(sig1_pk_all, pk_id_tmp)
bg_id_tmp = nbp_tmp>=fdr_thresh
frip_tmp = sum(sig_tmp[pk_id_tmp]) / sum(sig_tmp)
mean_pk = mean(sig_tmp[pk_id_tmp])
mean_bg = mean(sig_tmp[bg_id_tmp])
snr_tmp = mean_pk / mean_bg
pk_num_tmp = sum(pk_id_tmp)
ari_tmp = adjustedRandIndex(sig1_pk_id, pk_id_tmp)
ji_tmp = jaccard_index(sig1_pk_id, pk_id_tmp)
frip_common_tmp = frip_common(sig1, sig1_pk_id, sig_tmp, pk_id_tmp)
### append(info)
frip_all[dim(sig_matrix_p)[2]-1] = frip_tmp
snr_all[dim(sig_matrix_p)[2]-1] = snr_tmp
pk_num_all[dim(sig_matrix_p)[2]-1] = pk_num_tmp
ari_all[dim(sig_matrix_p)[2]-1] = ari_tmp
ji_all[dim(sig_matrix_p)[2]-1] = ji_tmp
frip_common_all = cbind(frip_common_all, frip_common_tmp)

### sig2z
nbp_tmp = sig_matrix_p[,dim(sig_matrix_p)[2]]
sig_tmp = sig_matrix[,dim(sig_matrix_p)[2]]
pk_id_tmp = nbp_tmp<fdr_thresh
sig1_pk_all = cbind(sig1_pk_all, pk_id_tmp)
bg_id_tmp = nbp_tmp>=fdr_thresh
frip_tmp = sum(sig_tmp[pk_id_tmp]) / sum(sig_tmp)
mean_pk = mean(sig_tmp[pk_id_tmp])
mean_bg = mean(sig_tmp[bg_id_tmp])
snr_tmp = mean_pk / mean_bg
pk_num_tmp = sum(pk_id_tmp)
ari_tmp = adjustedRandIndex(sig1_pk_id, pk_id_tmp)
ji_tmp = jaccard_index(sig1_pk_id, pk_id_tmp)
frip_common_tmp = frip_common(sig1, sig1_pk_id, sig_tmp, pk_id_tmp)
### append(info)
frip_all[dim(sig_matrix_p)[2]] = frip_tmp
snr_all[dim(sig_matrix_p)[2]] = snr_tmp
pk_num_all[dim(sig_matrix_p)[2]] = pk_num_tmp
ari_all[dim(sig_matrix_p)[2]] = ari_tmp
ji_all[dim(sig_matrix_p)[2]] = ji_tmp
frip_common_all = cbind(frip_common_all, frip_common_tmp)


### plot
ref_pk_id = sig_matrix_p[,dim(sig_matrix_p)[2]]<fdr_thresh
tar_pk_id = sig_matrix_p[,dim(sig_matrix_p)[2]-1]<fdr_thresh
cpk = (ref_pk_id * tar_pk_id) == 1
cbg = (ref_pk_id + tar_pk_id) == 0
ref_all_s = sig1_z[sample_id]
tar_all_s = sig2_z[sample_id]
cpk_s = cpk[sample_id]
cbg_s = cbg[sample_id]
ref_cpk_s = ref_all_s[cpk_s]
tar_cpk_s = tar_all_s[cpk_s]
ref_cbg_s = ref_all_s[cbg_s]
tar_cbg_s = tar_all_s[cbg_s]
cpk_mean = c(mean(sig1_z[cpk]), mean(sig2_z[cpk]))
cbg_mean = c(mean(sig1_z[cbg]), mean(sig2_z[cbg]))
print(cpk_mean)
print(cbg_mean)
all_mean = c(mean(sig1_z), mean(sig2_z))
print(all_mean)
png(paste(output_name, sig_matrix_colnames[dim(sig_matrix)[2]], '.scatter.png', sep=''))
plot_scatterplot_3parts(ref_all_s, tar_all_s, ref_cpk_s, tar_cpk_s, ref_cbg_s, tar_cbg_s, cpk_mean, cbg_mean, all_mean, -2, 7)
dev.off()



bed_file = as.data.frame(read.table(bed_file, header=F, sep='\t'))
all_pk_id = as.data.frame(sig1_pk_all*1)
colnames(all_pk_id) = sig_matrix_colnames
all_pk_id_sig = as.data.frame(cbind(all_pk_id, sig_matrix))
all_pk_bed = cbind(bed_file, all_pk_id_sig)
write.table(all_pk_bed, paste(output_name, '.pkid.txt', sep=''), quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')


info_matrix = cbind(frip_all, snr_all, pk_num_all, ari_all, ji_all, t(frip_common_all))
colnames(info_matrix) = c('frip', 'snr', 'pk_num', 'ari', 'ji', 'frip_cpk', 'frip_cpk_ref', 'frip_cpk_tar')
rownames(info_matrix) = sig_matrix_colnames

write.table(info_matrix, output_name, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')



