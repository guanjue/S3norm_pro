get_2r_nbp = function(sig, output_name){
	thesh = -1

	sig_non0 = sig[sig>=thesh]
	sig_mean = mean(sig_non0)
	sig_var = var(sig_non0)
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
	nb_pval = apply(as.matrix(sig_non0), MARGIN=1, function(x) pnbinom(x[1], sig_size, sig_prob, lower.tail=FALSE) )

	### 2nd round
	sig_non0_bg = sig_non0[nb_pval>=0.001]
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

	### write output
	if(!file.exists(output_name){
		write.table(nb_pval_2nd, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
	}
	return(nb_pval_2nd)
}


plot_MA_3parts = function(A_all, M_all, A_cpk, M_cpk, A_cbg, M_cbg, A_lim_lower, A_lim_upper){
	plot(A_all, M_all, pch=16, ylim=c(-M_ylim, M_ylim), xlim=c(A_lim_lower, A_lim_upper), col='dodgerblue', cex=0.6)
	points(A_cpk, M_cpk, pch=16, col='darkorange1', cex=0.6)
	points(A_cbg, M_cbg, pch=16, col='gray56', cex=0.6)
	points(mean(A_cpk), mean(M_cpk), pch=16, col='black', cex=1)
	points(mean(A_cbg), mean(M_cbg), pch=16, col='black', cex=1)
	lines(c(mean(A_cbg), mean(A_cpk)), c(mean(M_cbg), mean(M_cpk)), col='green', lty=2, lwd=3)
	abline(h=0, col='black', lwd=3)
}
