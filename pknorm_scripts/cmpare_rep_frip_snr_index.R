library(mclust)
library(Hmisc)

### get parameters
args = commandArgs(trailingOnly=TRUE)

sig1_raw = args[1]
sig2_raw = args[2]

sig1_TRnorm = args[3]
sig2_TRnorm = args[4]

sig1_MAnorm = args[5]
sig2_MAnorm = args[6]

sig1_QTnorm = args[7]
sig2_QTnorm = args[8]

sig1_PKnorm = args[9]
sig2_PKnorm = args[10]

output_name = args[11]
fdr_thresh = as.numeric(args[12])
method=args[13]
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


#sig1_raw = 'ER4.fisher_p.txt'
#sig2_raw = 'B_SPL.fisher_p.txt'
#sig2_TRnorm = 'B_SPL_TM_MAnorm_totalmean_MAnorm/B_SPL_TM_MAnorm.totalsig_norm.txt'
#sig2_MAnorm = 'B_SPL_TM_MAnorm_totalmean_MAnorm/B_SPL_TM_MAnorm.MAnorm.txt'
#sig2_QTnorm = 'B_SPL.fisher_p.txt.qtn.txt'
#sig2_PKnorm = 'B_SPL_txt.pknorm.txt'
#output_name = 'test'
#fdr_thresh = 0.05

### read data
sig1_r = scan(sig1_raw)
sig1_tr = scan(sig1_TRnorm)
sig1_ma = scan(sig1_MAnorm)
sig1_qt = scan(sig1_QTnorm)
sig1_pk = scan(sig1_PKnorm)

sig2_r = scan(sig2_raw)
sig2_tr = scan(sig2_TRnorm)
sig2_ma = scan(sig2_MAnorm)
sig2_qt = scan(sig2_QTnorm)
sig2_pk = scan(sig2_PKnorm)

sig_matrix = as.matrix(cbind(sig1_r, sig1_tr, sig1_ma, sig1_qt, sig1_pk, sig2_r, sig2_tr, sig2_ma, sig2_qt, sig2_pk))
sig_matrix_colnames = c('sig1_r', 'sig1_tr', 'sig1_ma', 'sig1_qt', 'sig1_pk', 'sig2_r', 'sig2_tr', 'sig2_ma', 'sig2_qt', 'sig2_pk')

### get nbp
if (method == 'nbp'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(nbp_2r(x, 0.001, paste(sig1_raw, '.nbp_2r.txt', sep='')), method='fdr'))
} else if (method == 'p'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(x, method='fdr'))	
} else if (method == 'neglog10p'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(10^(-x), method='fdr'))
} else if (method == 'zp'){
	sig_matrix_p = apply(sig_matrix, 2, function(x) p.adjust(zp(x), method='fdr'))
}

frip1_all = c()
snr1_all = c()
pk_num1_all = c()
frip2_all = c()
snr2_all = c()
pk_num2_all = c()
ari_all = c()
ji_all = c()

frip_common_all = c()

sig1_pk_id = sig_matrix_p[,1]<fdr_thresh
sig1 = sig_matrix[,1]

for ( i in c(1:(dim(sig_matrix)[2])/2)){
	print(sig_matrix_colnames[i])
	### get info
	nbp1_tmp = sig_matrix_p[,i]
	sig1_tmp = sig_matrix[,i]
	nbp2_tmp = sig_matrix_p[,i+dim(sig_matrix)[2])/2]
	sig2_tmp = sig_matrix[,i+dim(sig_matrix)[2])/2]
	pk1_id_tmp = nbp1_tmp<fdr_thresh
	bg1_id_tmp = nbp1_tmp>=fdr_thresh
	pk2_id_tmp = nbp2_tmp<fdr_thresh
	bg2_id_tmp = nbp2_tmp>=fdr_thresh
	frip1_tmp = sum(sig1_tmp[pk1_id_tmp]) / sum(sig1_tmp)
	snr1_tmp = mean(sig1_tmp[pk1_id_tmp]) / mean(sig1_tmp[bg1_id_tmp])
	pk_num1_tmp = sum(pk1_id_tmp)
	frip2_tmp = sum(sig2_tmp[pk2_id_tmp]) / sum(sig2_tmp)
	snr2_tmp = mean(sig2_tmp[pk2_id_tmp]) / mean(sig2_tmp[bg2_id_tmp])
	pk_num2_tmp = sum(pk2_id_tmp)
	ari_tmp = adjustedRandIndex(pk1_id_tmp, pk2_id_tmp)
	ji_tmp = jaccard_index(pk1_id_tmp, pk2_id_tmp)
	frip_common_tmp = frip_common(sig1_tmp, pk1_id_tmp, sig2_tmp, pk2_id_tmp)
	### append(info)
	frip1_all[i] = frip1_tmp
	snr1_all[i] = snr1_tmp
	pk_num1_all[i] = pk_num1_tmp
	frip2_all[i] = frip2_tmp
	snr2_all[i] = snr2_tmp
	pk_num2_all[i] = pk_num2_tmp
	ari_all[i] = ari_tmp
	ji_all[i] = ji_tmp
	frip_common_all = cbind(frip_common_all, frip_common_tmp)
}


info_matrix = cbind(frip1_all, snr1_all, pk_num1_all, frip2_all, snr2_all, pk_num2_all, ari_all, ji_all, t(frip_common_all))
colnames(info_matrix) = c('frip1', 'snr1', 'pk_num1', 'frip2', 'snr2', 'pk_num2', 'ari', 'ji', 'frip_cpk', 'frip_cpk_ref', 'frip_cpk_tar')
rownames(info_matrix) = sig_matrix_colnames

write.table(info_matrix, output_name, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')



