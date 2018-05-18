library(mclust)
library(Hmisc)

### get parameters
args = commandArgs(trailingOnly=TRUE)

sig1_raw = args[1]
sig2_raw = args[2]

sig1_QTnorm = args[3]
sig2_QTnorm = args[4]

sig2_PKnorm = args[5]
sig2_PKnorm_weight = args[6]
sig2_PKnorm_idx = args[7]

sig2_TRnorm = args[8]
sig2_MAnorm = args[9]


output_name = args[10]
fdr_thresh = as.numeric(args[11])

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
	write.table(nb_pval_2nd, output_name, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')
	### return 2nd round NB-p-value
	return(nb_pval_2nd)
}

################################################
plot_dif_col = function(sig_y, sig_x, common_pk, common_bg, min_sig, max_sig, output_name){
	pdf(output_name, width = 8, height = 8)
	plot(sig_x, sig_y, xlim=c(min_sig, max_sig), ylim=c(min_sig, max_sig), pch=16, cex=1, col = 'dodgerblue')
	points(sig_x[common_pk], sig_y[common_pk], col='darkorange1', pch=16, cex=1)
	points(sig_x[common_bg], sig_y[common_bg], col='gray28', pch=16, cex=1)
	points(mean(sig_x[common_pk]), mean(sig_y[common_pk]), col='black', pch=16, cex=2)
	points(mean(sig_x[common_bg]), mean(sig_y[common_bg]), col='black', pch=16, cex=2)
	points(mean(sig_x), mean(sig_y), col='red', pch=16, cex=2)
	abline(0,1,lwd=3,col='black')
	lines(c(mean(sig_x[common_bg]), mean(sig_x[common_pk])), c(mean(sig_y[common_bg]), mean(sig_y[common_pk])), col='royalblue1', lty=2, lwd=3)
	dev.off()
}

################################################
plot_dif_col_PKnorm = function(sig_y, sig_x, sig_x_weight, common_pk, common_bg, min_sig, max_sig, output_name){
	pdf(output_name, width = 8, height = 8)
	print('print PKnorm')
	print(head(cbind(sig_x[common_pk], sig_x_weight[common_pk])))
	plot(sig_x, sig_y, xlim=c(min_sig, max_sig), ylim=c(min_sig, max_sig), pch=16, cex=1, col = 'dodgerblue')
	points(sig_x[common_pk], sig_y[common_pk], col='darkorange1', pch=16, cex=1)
	points(sig_x[common_bg], sig_y[common_bg], col='gray28', pch=16, cex=1)
	points(wtd.mean(sig_x[common_pk], sig_x_weight[common_pk], normwt = "ignored"), wtd.mean(sig_y[common_pk], sig_x_weight[common_pk], normwt = "ignored"), col='black', pch=16, cex=2)
	points(wtd.mean(sig_x[common_bg], sig_x_weight[common_bg], normwt = "ignored"), wtd.mean(sig_y[common_bg], sig_x_weight[common_bg], normwt = "ignored"), col='black', pch=16, cex=2)
	points(mean(sig_x), mean(sig_y), col='red', pch=16, cex=2)
	abline(0,1,lwd=3,col='black')
	lines(c(wtd.mean(sig_x[common_bg], sig_x_weight[common_bg], normwt = "ignored"), wtd.mean(sig_x[common_pk], sig_x_weight[common_pk], normwt = "ignored")), c(wtd.mean(sig_y[common_bg], sig_x_weight[common_bg], normwt = "ignored"), wtd.mean(sig_y[common_pk], sig_x_weight[common_pk], normwt = "ignored")), col='royalblue1', lty=2, lwd=3)
	dev.off()
}

################################################
get_FRiP_pknum_jaccard_index = function(sig_y, sig_x, sig_y_nbp, sig_x_nbp){
	### FRiP
	sig_x_FRiP = sum(sig_x[sig_x_nbp]) / sum(sig_x)
	sig_y_FRiP = sum(sig_y[sig_y_nbp]) / sum(sig_y)
	### pk & bg mean
	sig_x_pk_mean = mean(sig_x[sig_x_nbp])
	sig_y_pk_mean = mean(sig_y[sig_y_nbp])
	sig_x_bg_mean = mean(sig_x[sig_x_nbp==0])
	sig_y_bg_mean = mean(sig_y[sig_y_nbp==0])
	sig_x_total_mean = mean(sig_x)
	sig_y_total_mean = mean(sig_y)

	### pk_num
	sig_x_pknum = sum(sig_x_nbp)
	sig_y_pknum = sum(sig_y_nbp)
	### jaccard_index
	overlap = (sig_x_nbp * sig_y_nbp)==1
	pk_overlap_num = sum(overlap)
	union = (sig_x_nbp + sig_y_nbp)!=0
	pk_union_num = sum(union)
	jaccard_index = pk_overlap_num / pk_union_num

	FRiP_JI_all = (sum(sig_x[overlap])+sum(sig_y[overlap])) / (sum(sig_x)+sum(sig_y))
	FRiP_JI_x = sum(sig_x[overlap]) / sum(sig_x)
	FRiP_JI_y = sum(sig_y[overlap]) / sum(sig_y)

	uniq_x = (sig_x_nbp * (!sig_y_nbp))==1
	uniq_y = (sig_y_nbp * (!sig_x_nbp))==1
	uniq_xy = (uniq_x + uniq_y) ==1
	FRuP_JI_all = (sum(sig_x[uniq_x])+sum(sig_y[uniq_y])) / (sum(sig_x)+sum(sig_y))
	FRuP_JI_x = sum(sig_x[uniq_x]) / sum(sig_x)
	FRuP_JI_y = sum(sig_y[uniq_y]) / sum(sig_y)


	ari = adjustedRandIndex(sig_x_nbp, sig_y_nbp)

	return(c(sig_x_FRiP, sig_y_FRiP, sig_x_pknum, sig_y_pknum, jaccard_index, ari, pk_overlap_num, pk_union_num, sig_x_pk_mean, sig_y_pk_mean, sig_x_bg_mean, sig_y_bg_mean, sig_x_total_mean, sig_y_total_mean, FRiP_JI_all, FRiP_JI_x, FRiP_JI_y, FRuP_JI_all, FRuP_JI_x, FRuP_JI_y))
}

################################################
mvfile2folder = function(from, to) {
    todir = dirname(to)
    if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
    file.rename(from = from,  to = to)
}

################################################
plot_5 = function(d1_raw, d2_raw, d1_QTnorm, d2_PKnorm, d2_PKnorm_weight, d2_PKnorm_idx, d2_TRnorm, d2_MAnorm, d2_QTnorm, d1_raw_nbp, d2_raw_nbp, d1_QTnorm_nbp, d2_PKnorm_nbp, d2_TRnorm_nbp, d2_MAnorm_nbp, d2_QTnorm_nbp, output_name){
	### convert to log2 scale
	d1_raw_log2 = log2(d1_raw+0.1)
	d2_raw_log2 = log2(d2_raw+0.1)
	d1_QTnorm_log2 = log2(d1_QTnorm+0.1)
	d2_PKnorm_log2 = log2(d2_PKnorm+0.1)
	d2_TRnorm_log2 = log2(d2_TRnorm+0.1)
	d2_MAnorm_log2 = log2(d2_MAnorm+0.1)
	d2_QTnorm_log2 = log2(d2_QTnorm+0.1)
	### read PKnorm weight
	d2_PKnorm_weight = d2_PKnorm_weight

	### sampling for plotting
	set.seed(2018)
	sample_id = d2_PKnorm_idx
	### sample signals
	print('sample signals!!!')
	d1_raw_log2_s = d1_raw_log2[sample_id]
	d2_raw_log2_s = d2_raw_log2[sample_id]
	d1_QTnorm_log2_s = d1_QTnorm_log2[sample_id]
	d2_PKnorm_log2_s = d2_PKnorm_log2[sample_id]
	d2_TRnorm_log2_s = d2_TRnorm_log2[sample_id]
	d2_MAnorm_log2_s = d2_MAnorm_log2[sample_id]
	d2_QTnorm_log2_s = d2_QTnorm_log2[sample_id]

	d2_PKnorm_weight_s = d2_PKnorm_weight

	### sample NB-p
	d1_raw_nbp_s = d1_raw_nbp[sample_id]
	d2_raw_nbp_s = d2_raw_nbp[sample_id]
	d1_QTnorm_nbp_s = d1_QTnorm_nbp[sample_id]
	d2_PKnorm_nbp_s = d2_PKnorm_nbp[sample_id]
	d2_TRnorm_nbp_s = d2_TRnorm_nbp[sample_id]
	d2_MAnorm_nbp_s = d2_MAnorm_nbp[sample_id]
	d2_QTnorm_nbp_s = d2_QTnorm_nbp[sample_id]

	### get max and min
	all_matrix_sample = cbind(d1_raw_log2_s, d2_raw_log2_s, d1_QTnorm_log2_s, d2_PKnorm_log2_s, d2_TRnorm_log2_s, d2_MAnorm_log2_s, d2_QTnorm_log2_s)
	min_sig = min(all_matrix_sample)
	max_sig = max(all_matrix_sample)

	### plot figure
	d12_raw_nbp_common_pk_s = (d1_raw_nbp_s * d2_raw_nbp_s) == 1
	d12_raw_nbp_common_bg_s = (d1_raw_nbp_s + d2_raw_nbp_s) == 0
	plot_dif_col(d1_raw_log2_s, d2_raw_log2_s, d12_raw_nbp_common_pk_s, d12_raw_nbp_common_bg_s, min_sig, max_sig, paste(output_name, '.raw.pdf', sep=''))

	d12_QTnorm_nbp_common_pk_s = (d1_QTnorm_nbp_s * d2_QTnorm_nbp_s) == 1
	d12_QTnorm_nbp_common_bg_s = (d1_QTnorm_nbp_s + d2_QTnorm_nbp_s) == 0
	plot_dif_col(d1_QTnorm_log2_s, d2_QTnorm_log2_s, d12_QTnorm_nbp_common_pk_s, d12_QTnorm_nbp_common_bg_s, min_sig, max_sig, paste(output_name, '.QTnorm.pdf', sep=''))

	d12_PKnorm_nbp_common_pk_s = (d1_raw_nbp_s * d2_PKnorm_nbp_s) == 1
	d12_PKnorm_nbp_common_bg_s = (d1_raw_nbp_s + d2_PKnorm_nbp_s) == 0
	print(head(d12_PKnorm_nbp_common_pk_s))
	print(length(d2_PKnorm_log2_s))
	print(length(d2_PKnorm_weight_s))
	#plot_dif_col_PKnorm(d1_raw_log2_s, d2_PKnorm_log2_s, d2_PKnorm_weight_s, d12_PKnorm_nbp_common_pk_s, d12_PKnorm_nbp_common_bg_s, min_sig, max_sig, paste(output_name, '.PKnorm.pdf', sep=''))
	plot_dif_col(d1_raw_log2_s, d2_PKnorm_log2_s, d12_PKnorm_nbp_common_pk_s, d12_PKnorm_nbp_common_bg_s, min_sig, max_sig, paste(output_name, '.PKnorm.pdf', sep=''))

	d12_TRnorm_nbp_common_pk_s = (d1_raw_nbp_s * d2_TRnorm_nbp_s) == 1
	d12_TRnorm_nbp_common_bg_s = (d1_raw_nbp_s + d2_TRnorm_nbp_s) == 0
	plot_dif_col(d1_raw_log2_s, d2_TRnorm_log2_s, d12_TRnorm_nbp_common_pk_s, d12_TRnorm_nbp_common_bg_s, min_sig, max_sig, paste(output_name, '.TRnorm.pdf', sep=''))

	d12_MAnorm_nbp_common_pk_s = (d1_raw_nbp_s * d2_MAnorm_nbp_s) == 1
	d12_MAnorm_nbp_common_bg_s = (d1_raw_nbp_s + d2_MAnorm_nbp_s) == 0
	plot_dif_col(d1_raw_log2_s, d2_MAnorm_log2_s, d12_MAnorm_nbp_common_pk_s, d12_MAnorm_nbp_common_bg_s, min_sig, max_sig, paste(output_name, '.MAnorm.pdf', sep=''))


	dir.create(paste(output_name, '_difnorm_compare', sep=''), showWarnings = FALSE)
	mvfile2folder(from=paste(output_name, '.raw.pdf', sep=''), to=paste(output_name, '_difnorm_compare/', output_name, '.raw.pdf', sep=''))
	mvfile2folder(from=paste(output_name, '.QTnorm.pdf', sep=''), to=paste(output_name, '_difnorm_compare/', output_name, '.QTnorm.pdf', sep=''))
	mvfile2folder(from=paste(output_name, '.PKnorm.pdf', sep=''), to=paste(output_name, '_difnorm_compare/', output_name, '.PKnorm.pdf', sep=''))
	mvfile2folder(from=paste(output_name, '.TRnorm.pdf', sep=''), to=paste(output_name, '_difnorm_compare/', output_name, '.TRnorm.pdf', sep=''))
	mvfile2folder(from=paste(output_name, '.MAnorm.pdf', sep=''), to=paste(output_name, '_difnorm_compare/', output_name, '.MAnorm.pdf', sep=''))
}
################################################


### read raw signal
print('read raw files')
d1_raw = scan(sig1_raw)
d1_raw_nbp = p.adjust(nbp_2r(d1_raw, 0.001, paste(sig1_raw, '.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh
d2_raw = scan(sig2_raw)
d2_raw_nbp = p.adjust(nbp_2r(d2_raw, 0.001, paste(sig2_raw, '.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh

### read Quantile normalized signal
print('read QTnorm files')
d1_QTnorm = scan(sig1_QTnorm)
d1_QTnorm_nbp = p.adjust(nbp_2r(d1_QTnorm, 0.001, paste(sig1_raw, '.QTnorm.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh
d2_QTnorm = scan(sig2_QTnorm)
d2_QTnorm_nbp = p.adjust(nbp_2r(d2_QTnorm, 0.001, paste(sig2_raw, '.QTnorm.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh


### read PKnorm, total mean normalized, MAnorm
print('read PKnorm files')
d2_PKnorm = scan(sig2_PKnorm)
d2_PKnorm_nbp = p.adjust(nbp_2r(d2_PKnorm, 0.001, paste(sig2_raw, '.PKnorm.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh
d2_PKnorm_weight = scan(sig2_PKnorm_weight)
d2_PKnorm_idx = scan(sig2_PKnorm_idx)
print('read TRnorm files')
d2_TRnorm = scan(sig2_TRnorm)
d2_TRnorm_nbp = p.adjust(nbp_2r(d2_TRnorm, 0.001, paste(sig2_raw, '.TRnorm.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh
print('read MAnorm files')
d2_MAnorm = scan(sig2_MAnorm)
d2_MAnorm_nbp = p.adjust(nbp_2r(d2_MAnorm, 0.001, paste(sig2_raw, '.MAnorm.nbp_2r.txt', sep='')), method='fdr') < fdr_thresh

print('plot 5 files')
plot_5(d1_raw, d2_raw, d1_QTnorm, d2_PKnorm, d2_PKnorm_weight, d2_PKnorm_idx, d2_TRnorm, d2_MAnorm, d2_QTnorm, d1_raw_nbp, d2_raw_nbp, d1_QTnorm_nbp, d2_PKnorm_nbp, d2_TRnorm_nbp, d2_MAnorm_nbp, d2_QTnorm_nbp, output_name)

print('get numeric summary')
FRiP_pknum_JI_raw = get_FRiP_pknum_jaccard_index(d1_raw, d2_raw, d1_raw_nbp, d2_raw_nbp)
FRiP_pknum_JI_QTnorm = get_FRiP_pknum_jaccard_index(d1_QTnorm, d2_QTnorm, d1_QTnorm_nbp, d2_QTnorm_nbp)
FRiP_pknum_JI_PKnorm = get_FRiP_pknum_jaccard_index(d1_raw, d2_PKnorm, d1_raw_nbp, d2_PKnorm_nbp)
FRiP_pknum_JI_TRnorm = get_FRiP_pknum_jaccard_index(d1_raw, d2_TRnorm, d1_raw_nbp, d2_TRnorm_nbp)
FRiP_pknum_JI_MAnorm = get_FRiP_pknum_jaccard_index(d1_raw, d2_MAnorm, d1_raw_nbp, d2_MAnorm_nbp)

FRiP_pknum_JI_matrix = as.matrix(rbind(FRiP_pknum_JI_raw, FRiP_pknum_JI_QTnorm, FRiP_pknum_JI_PKnorm, FRiP_pknum_JI_TRnorm, FRiP_pknum_JI_MAnorm))
rownames(FRiP_pknum_JI_matrix) = c('raw', 'QTnorm', 'PKnorm', 'TRnorm', 'MAnorm')
colnames(FRiP_pknum_JI_matrix) = c('sig_x_FRiP', 'sig_y_FRiP', 'sig_x_pknum', 'sig_y_pknum', 'jaccard_index', 'adjusted_rand_index', 'pk_overlap_num', 'pk_union_num', 'sig_x_pk_mean', 'sig_y_pk_mean', 'sig_x_bg_mean', 'sig_y_bg_mean', 'sig_x_total_mean', 'sig_y_total_mean', 'FRiP_JI_all', 'FRiP_JI_x', 'FRiP_JI_y', 'FRuP_JI_all', 'FRuP_JI_x', 'FRuP_JI_y')

write.table(FRiP_pknum_JI_matrix, paste(output_name, '.summary_matrix.txt', sep=''), quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')
mvfile2folder(from=paste(output_name, '.summary_matrix.txt', sep=''), to=paste(output_name, '_difnorm_compare/', output_name, '.summary_matrix.txt', sep=''))

mvfile2folder(from=paste(sig1_raw, '.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig1_raw, '.nbp_2r.txt', sep=''))
mvfile2folder(from=paste(sig2_raw, '.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig2_raw, '.nbp_2r.txt', sep=''))
mvfile2folder(from=paste(sig1_raw, '.QTnorm.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig1_raw, '.QTnorm.nbp_2r.txt', sep=''))
mvfile2folder(from=paste(sig2_raw, '.QTnorm.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig2_raw, '.QTnorm.nbp_2r.txt', sep=''))
mvfile2folder(from=paste(sig2_raw, '.PKnorm.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig2_raw, '.PKnorm.nbp_2r.txt', sep=''))
mvfile2folder(from=paste(sig2_raw, '.TRnorm.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig2_raw, '.TRnorm.nbp_2r.txt', sep=''))
mvfile2folder(from=paste(sig2_raw, '.MAnorm.nbp_2r.txt', sep=''), to=paste(output_name, '_difnorm_compare/', sig2_raw, '.MAnorm.nbp_2r.txt', sep=''))


