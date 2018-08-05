#library(robustreg)
library(MASS)
library(affy)

get_nbp = function(sig){
	thesh = 0
	sig_non0 = sig[sig>thesh]
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
	return(nb_pval_2nd)
}


args = commandArgs(trailingOnly=TRUE)
input_sig1 = args[1]
input_sig2 = args[2]
output = args[3]

small_num = 0.01
random_sample_num = 1000000
upperlim = 1000
lowerlim = 0
ct = unlist(strsplit(input_sig2, "[.]"))[1]

sig1 = scan(input_sig1)
sig2 = scan(input_sig2)
totalmean_sf = sum(sig1) / sum(sig2)
sig3 = (sig2) * totalmean_sf #- small_num
sig4 = scale(sig2) #- small_num


#sig1_binary = get_nbp(sig1) <= 0.001
#sig2_binary = get_nbp(sig2) <= 0.001

sig1_binary = 10^(-sig1) <= 0.001
sig2_binary = 10^(-sig2) <= 0.001

#sig1_binary = p.adjust(10^(-sig1),'fdr') < 0.05
if (sum(sig1_binary)<as.integer(dim(sig1)[1]/100)){
	pk1_lim = sort(sig2)[length(sig1)-as.integer(dim(sig1)[1]/100)]
	sig1_binary = sig1 >= pk1_lim
}
#sig2_binary = p.adjust(10^(-sig2),'fdr') < 0.05
if (sum(sig2_binary)<as.integer(dim(sig1)[1]/100)){
	pk2_lim = sort(sig2)[length(sig2)-as.integer(dim(sig1)[1]/100)]
	sig2_binary = sig2 >= pk2_lim
}

peak_binary_pk = as.logical(sig1_binary * sig2_binary) 
peak_binary = peak_binary_pk & (sig1 != sig1[1]) & (sig2 != sig2[1])

bg_binary_bg = as.logical((sig1_binary + sig2_binary)==0)
bg_binary = bg_binary_bg & (sig1 != sig1[1]) & (sig2 != sig2[1]) 

common_peak_count_read1 = sig1[peak_binary]
common_peak_count_read2 = sig2[peak_binary]


M=log2((common_peak_count_read2+small_num)/(common_peak_count_read1+small_num))
A=0.5*log2((common_peak_count_read2+small_num)*(common_peak_count_read1+small_num))
M = as.matrix(M)
A = as.matrix(A)

linear=lm(M~A)$coefficients
#b=lm(M~A)$coefficients
#b=robustRegBS(M,A,beta=linear)
b=rlm(M~A)$coefficients

cat("M = b[1] + b[2] * A\n")

log2_allregion_count_read1 = log2(sig1 + small_num)
log2_allregion_count_read2 = log2(sig2 + small_num)
log2_allregion_count_read2_rescaled = (2-b[2])*log2_allregion_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
sig2_rescaled = 2^log2_allregion_count_read2_rescaled - small_num

#write.table(info, paste(output, '.MA.norm.info.txt', sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(sig2_rescaled, paste(output,".manorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(sig3, paste(output,".trnorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(sig4, paste(output,".znorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)





