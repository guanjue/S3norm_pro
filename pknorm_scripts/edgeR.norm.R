library(edgeR)
library(seewave)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(LSD)

d = read.table('rnaHtseqCountsall_replicate_merge.pcsorted.txt', header=FALSE)
d_sig = d[,-c(1,3)]
ct_list = c('CFU_E_ad', 'CMP', 'ERY_ad', 'GMP', 'MK_imm_ad', 'LSK_BM', 'MEP', 'MONO_BM', 'NEU', 'ER4', 'G1E')
colnames(d_sig) = ct_list
dgList <- DGEList(counts=d_sig, genes=d[,1])
### filter genes
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
#dgList <- dgList[keep,]
#summary(cpm(dgList)) #compare this to the original summary

dgList <- calcNormFactors(dgList, method="TMM")
countsPerMillion_filter <- cpm(dgList)

regions = read.table('gencode.vM4.annotation.pc.sorted.bed', header=FALSE)
regions_length = regions[,3]-regions[,2]

library(LSD)
heatscatter(countsPerMillion_filter[,1],countsPerMillion_filter[,6], log='xy')
abline(0,1,col='red')

heatscatter(data0[,2],d_sig[,1], log='xy')
abline(0,1,col='red')
heatscatter(data0[,3],d[,3], log='xy')
abline(0,1,col='red')
heatscatter(data0[,4],d[,4], log='xy')
abline(0,1,col='red')


used_id_rna_tpm = keep

rna_tpm = (countsPerMillion_filter / regions_length)[keep,]

methods = c('rcnorm', 'poisnorm', 'rcznorm', 'raw', 'trnorm', 'qtnorm', 'manorm', 'pknorm')
#methods = c('rcnorm', 'poisnorm', 'rcznorm', 'poisnorm', 'raw', 'trnorm', 'qtnorm', 'manorm', 'pknorm')

#methods = c('rcnorm', 'rcznorm', 'raw', 'znorm', 'trnorm', 'qtnorm', 'manorm', 'pknorm')

small_num = 0.1


cor_0_matrix = c()
cor_0_shuffle_matrix = c()
paired_t_statistic_vec = c()
kl_dist_vec = c()

set.seed(2018)
for (m in methods){
	print(m)
	set.seed(1)
	bw_used = 0.1
	cor_method = 'spearman'
	#small_num = 0.1
	#shuffle_id = sample(dim(rna_tpm)[1],dim(rna_tpm)[1])
	if (m!='rcznorm'){
		d_raw = (as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))+small_num)
		d_raw = ((d_raw) - mean((d_raw))) / sd((d_raw))
	} else {
		d_raw = (as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F)))
	}
	print(summary(d_raw))
	d_raw = d_raw[keep,]
	
	#d_raw = log2(d_raw+small_num)
	sig_mat = cbind(rna_tpm, d_raw)
	cor_0 = apply(sig_mat, 1, function(x) cor((x[1:11]),(x[12:22]), method=cor_method))
	shuffle_id = sample(dim(rna_tpm)[2],dim(rna_tpm)[2])
	#d_raw_shuffle = d_raw[,shuffle_id]
	d_raw_shuffle = t(apply(d_raw, 1, function(x) x[sample(dim(rna_tpm)[2],dim(rna_tpm)[2])]))
	cor_0_shuffle = apply(cbind(rna_tpm, d_raw_shuffle), 1, function(x) cor((x[1:11]),(x[12:22]), method=cor_method))
	#d_raw_bg = read.table(paste('tss_h3k4me3.pcsorted.', m, '.1000kb.txt', sep=''), header=F)
	#cor_0_bg = apply(cbind(rna_tpm, d_raw_bg), 1, function(x) cor((x[1:11]),(x[12:22]), method=cor_method))
	#kl_dist_bg = kl.dist(density(cor_0[!is.na(cor_0)], bw=bw_used)$y, density(cor_0_bg[!is.na(cor_0_bg)], bw=bw_used)$y)$D2
	kl_dist_shuffle = kl.dist(density(cor_0[!is.na(cor_0)], bw=bw_used)$y, density(cor_0_shuffle[!is.na(cor_0_shuffle)], bw=bw_used)$y)$D2
	paired_t = t.test(cor_0, cor_0_shuffle, paired=TRUE, alternative = 'greater')
	paired_t_statistic = paired_t$statistic	
	print(summary(cor_0[!is.na(cor_0)]))
	png(paste('tss_h3k4me3.pcsorted.', m, '.png', sep=''))
	plot(density(cor_0[!is.na(cor_0)], bw=bw_used), col='green', main=paste('Paired-t-test-statistic = ', toString(round(paired_t_statistic, digits=3)), '; ', 'KL-dist = ', toString(round(kl_dist_shuffle, digits=3)), sep=''), ylim=c(0,1.5))
	print(mean(cor_0[!is.na(cor_0)]))
	print(median(cor_0[!is.na(cor_0)]))
	#plot(density(cor_0[!is.na(cor_0)], bw=bw_used), col='green', main=paste('Paired-t-test-statistic = ', toString(round(paired_t_statistic, digits=3)), sep=''), ylim=c(0,1.2))
	#lines(density(cor_0_bg[!is.na(cor_0_bg)], bw=bw_used), col='black')
	lines(density(cor_0_shuffle[!is.na(cor_0_shuffle)], bw=bw_used), col='blue')
	dev.off()
	cor_0_matrix = cbind(cor_0_matrix, cor_0[!is.na(cor_0)])
	cor_0_shuffle_matrix = cbind(cor_0_shuffle_matrix, cor_0_shuffle[!is.na(cor_0_shuffle)])
	paired_t_statistic_vec = cbind(paired_t_statistic_vec, paired_t_statistic)
	kl_dist_vec = cbind(kl_dist_vec, kl_dist_shuffle)
}





hist(rna_tpm[cor_0>0.5,])

hist(rna_tpm)


pdf('cor_0_KL-dist.pdf', width=20, height=9)
#barplot(kl_dist_vec[,c(1,2,3,4,5,6,7,8,9,10,11,12)], main="KL-dist", xlab="Methods", names.arg=c('RC', 'RC-Z', 'RC-TSnorm', 'RC-QTnorm', 'RC-MAnorm', 'POISP', 'NBP', 'NBP-Z', 'NBP-TSnorm', 'NBP-QTnorm', 'NBP-MAnorm', 'NBP-PKnorm'))
barplot(kl_dist_vec[,c(1,2,3,4,5,6,7,8)], main="KL-dist", xlab="Methods", names.arg=c('RC', 'POISP', 'RC-Z', 'NBP', 'NBP-TSnorm', 'NBP-QTnorm', 'NBP-MAnorm', 'NBP-PKnorm'))

#barplot(kl_dist_vec[,c(1,2,3,4,5,6,7,8)], main="KL-dist", xlab="Methods", names.arg=c('RC', 'RC-Z', 'NBP', 'NBP-Z', 'NBP-TSnorm', 'NBP-QTnorm', 'NBP-MAnorm', 'NBP-PKnorm'))
#box()
dev.off()

pdf('cor_0_paired_t_statistic_vec.pdf', width=20, height=9)
#barplot(paired_t_statistic_vec[,c(1,2,3,4,5,6,7,8,9,10,11,12)], main="KL-dist", xlab="Methods", names.arg=c('RC', 'RC-Z', 'RC-TSnorm', 'RC-QTnorm', 'RC-MAnorm', 'POISP', 'NBP', 'NBP-Z', 'NBP-TSnorm', 'NBP-QTnorm', 'NBP-MAnorm', 'NBP-PKnorm'))
barplot(paired_t_statistic_vec[,c(1,2,3,4,5,6,7,8)], main="paired_t_statistic_vec", xlab="Methods", names.arg=c('RC', 'POISP', 'RC-Z', 'NBP', 'NBP-TSnorm', 'NBP-QTnorm', 'NBP-MAnorm', 'NBP-PKnorm'))

#barplot(paired_t_statistic_vec[,c(1,2,3,4,5,6,7,8)], main="KL-dist", xlab="Methods", names.arg=c('RC', 'RC-Z', 'NBP', 'NBP-Z', 'NBP-TSnorm', 'NBP-QTnorm', 'NBP-MAnorm', 'NBP-PKnorm'))
#box()
dev.off()
