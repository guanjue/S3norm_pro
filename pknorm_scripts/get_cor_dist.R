library(seewave)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

data0 =read.table('rna_rpk.pcsorted.txt', header=F)
#data01 =read.table('rna_rpk.pcsorted.qtnorm.txt', header=F)

#data0 =read.table('rna_rpk.pcsorted.all12.txt', header=F)

total_rpk = colSums(data0[,-1])

rna_tpm = t(apply(data0[,-1], 1, function(x) x/total_rpk*10000000))

rna_tpm_max = apply(rna_tpm, 1, max)
rna_tpm_mean = apply(rna_tpm, 1, mean)
tpm_lim=5
sd_mean_lim = 1
non0_count_lim = 4
pdf('hist_rna_tpm_max.pdf')
hist(log2(rna_tpm_max), breaks=50)
abline(v=tpm_lim, col='red', lwd=1.5, lty=3)
dev.off()

rna_tpm_binary = (rna_tpm>tpm_lim)*1
rna_tpm_binary_count = colSums(rna_tpm_binary)

ct_list = c('CFU_E_ad', 'CMP', 'ERY_ad', 'GMP', 'MK_imm_ad', 'LSK_BM', 'MEP', 'MONO_BM', 'NEU', 'ER4', 'G1E')
ct_new_order = c(6,2,7,11,10,1,3,5,4,8,7)
ct_list_reorder = ct_list[ct_new_order]

png('rna_pk_count.png')
plot(c(1:length(rna_tpm_binary_count)), as.numeric(rna_tpm_binary_count)[ct_new_order], xaxt = 'n')
axis(1, at=c(1:length(rna_tpm_binary_count)), labels = FALSE)
text(c(1:length(rna_tpm_binary_count)), par("usr")[3], labels = ct_list_reorder, srt = 45, pos = 1, xpd = TRUE)
dev.off()




tpm_lim_mean = 3
pdf('hist_rna_tpm_mean.pdf')
hist(log2(rna_tpm_mean), breaks=50)
abline(v=tpm_lim_mean, col='red', lwd=1.5, lty=3)
dev.off()


rna_tpm_min = apply(rna_tpm, 1, min)
rna_tpm_min = apply(rna_tpm, 1, min)

rna_tpm_non0_num = apply(rna_tpm, 1, function(x) sum(x!=0))
rna_tpm_cv = apply(rna_tpm, 1, function(x) (sd(x+1))/mean(x+1))

#used_id_rna_tpm = ( (rna_tpm_max>=tpm_lim) * (rna_tpm_non0_num>=non0_count_lim) * (rna_tpm_cv>sd_mean_lim)) >0
#used_id_rna_tpm = ( (log2(rna_tpm_mean+0.01)>tpm_lim_mean) * (log2(rna_tpm_max+0.01)>tpm_lim) ) >0
#used_id_rna_tpm = ( (rna_tpm_cv>sd_mean_lim) * (log2(rna_tpm_max+0.01)>tpm_lim) ) >0
used_id_rna_tpm = ( (log2(rna_tpm_max+0.1)>tpm_lim)  ) >0
#used_id_rna_tpm = ( (rna_tpm_max>-1000)  ) >0
print(sum(used_id_rna_tpm))

#pdf('hist_tpm_cv.pdf', width=1200, height=400)
pdf('hist_tpm_cv.pdf')
#par(mfrow=c(1,3))
#hist(log2(rna_tpm+0.01), breaks=50)
#abline(v=tpm_lim, col='red', lwd=1.5, lty=3)
#box()
#hist(rna_tpm_non0_num, breaks=50)
#abline(v=non0_count_lim, col='red', lwd=1.5, lty=3)
#box()
hist(rna_tpm_cv, breaks=50)
abline(v=sd_mean_lim, col='red', lwd=1.5, lty=3)
box()
dev.off()

###### get cross cell type correlation
#methods = c('rcnorm', 'rcznorm', 'rctrnorm', 'rcqtnorm', 'rcmanorm', 'poisnorm', 'raw', 'znorm', 'trnorm', 'qtnorm', 'manorm', 'pknorm')
methods = c('rcnorm', 'poisnorm', 'rcznorm', 'raw', 'trnorm', 'qtnorm', 'manorm', 'pknorm')
#methods = c('rcnorm', 'poisnorm', 'rcznorm', 'poisnorm', 'raw', 'trnorm', 'qtnorm', 'manorm', 'pknorm')

#methods = c('rcnorm', 'rcznorm', 'raw', 'znorm', 'trnorm', 'qtnorm', 'manorm', 'pknorm')

small_num = 0.1


rna_tpm = as.matrix(log10(rna_tpm+small_num))
rna_tpm = (rna_tpm-mean(rna_tpm))/sd(rna_tpm)
#rna_tpm = scale(rna_tpm)
#rna_tpm = log2(rna_tpm+0.01)
print(summary(rna_tpm))



cor_0_matrix = c()
cor_0_shuffle_matrix = c()
paired_t_statistic_vec = c()
kl_dist_vec = c()

set.seed(2018)
for (m in methods){
	print(m)
	set.seed(1)
	bw_used = 0.1
	cor_method = 'pearson'
	#small_num = 0.1
	#shuffle_id = sample(dim(rna_tpm)[1],dim(rna_tpm)[1])
	if (m!='rcznorm'){
		d_raw = (as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))+small_num)
		d_raw = ((d_raw) - mean((d_raw))) / sd((d_raw))
	} else {
		d_raw = (as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F)))
	}
	print(summary(d_raw))
	
	#d_raw = log2(d_raw+small_num)
	sig_mat = cbind(rna_tpm, d_raw)[used_id_rna_tpm,]
	cor_0 = apply(sig_mat, 1, function(x) cor((x[1:11]),(x[12:22]), method=cor_method))
	shuffle_id = sample(dim(rna_tpm)[2],dim(rna_tpm)[2])
	#d_raw_shuffle = d_raw[,shuffle_id]
	d_raw_shuffle = t(apply(d_raw, 1, function(x) x[sample(dim(rna_tpm)[2],dim(rna_tpm)[2])]))
	cor_0_shuffle = apply(cbind(rna_tpm, d_raw_shuffle)[used_id_rna_tpm,], 1, function(x) cor((x[1:11]),(x[12:22]), method=cor_method))
	d_raw_bg = read.table(paste('tss_h3k4me3.pcsorted.', m, '.1000kb.txt', sep=''), header=F)
	cor_0_bg = apply(cbind(rna_tpm, d_raw_bg)[used_id_rna_tpm,], 1, function(x) cor((x[1:11]),(x[12:22]), method=cor_method))
	kl_dist_bg = kl.dist(density(cor_0[!is.na(cor_0)], bw=bw_used)$y, density(cor_0_bg[!is.na(cor_0_bg)], bw=bw_used)$y)$D2
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



library(LSD)
png('check.png')
heatscatter(rna_tpm[used_id_rna_tpm,1],rna_tpm[used_id_rna_tpm,6])
print(mean(rna_tpm[used_id_rna_tpm,1]))
print(mean(rna_tpm[used_id_rna_tpm,6]))
print(median(rna_tpm[used_id_rna_tpm,1]))
print(median(rna_tpm[used_id_rna_tpm,6]))
abline(0,1,col='red')
dev.off()

png('check1.png')
heatscatter(log10(data0[,-1]+0.1)[used_id_rna_tpm,1],log10(data0[,-1]+0.1)[used_id_rna_tpm,6])
abline(0,1,col='red')
dev.off()


heatscatter(d_raw[used_id_rna_tpm,1],d_raw[used_id_rna_tpm,6])
abline(0,1,col='red')













png('check_raw.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'raw', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
#heatscatter(rna_tpm[rna_tpm[,colnum]!=0,colnum]-rna_tpm[rna_tpm[,colnum]!=0,colnum+1], d_raw[rna_tpm[,colnum]!=0,colnum]-d_raw[rna_tpm[,colnum]!=0,colnum+1], log='', main=toString(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson')))
heatscatter(rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1], d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1], log='', main=toString(sum(((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])-(d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1]))^2)/sum((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])^2)) )

abline(0,1,col='red')
print(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson'))
dev.off()

png('check_qt.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'qtnorm', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
#heatscatter(rna_tpm[rna_tpm[,colnum]!=0,colnum]-rna_tpm[rna_tpm[,colnum]!=0,colnum+1], d_raw[rna_tpm[,colnum]!=0,colnum]-d_raw[rna_tpm[,colnum]!=0,colnum+1], log='', main=toString(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson')))
heatscatter(rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1], d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1], log='', main=toString(sum(((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])-(d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1]))^2)/sum((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])^2)) )

abline(0,1,col='red')
print(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson'))
dev.off()


png('check_ma.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'manorm', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
#heatscatter(rna_tpm[rna_tpm[,colnum]!=0,colnum]-rna_tpm[rna_tpm[,colnum]!=0,colnum+1], d_raw[rna_tpm[,colnum]!=0,colnum]-d_raw[rna_tpm[,colnum]!=0,colnum+1], log='', main=toString(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson')))
heatscatter(rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1], d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1], log='', main=toString(sum(((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])-(d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1]))^2)/sum((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])^2)) )

abline(0,1,col='red')
print(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson'))
dev.off()

png('check_ts.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'trnorm', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
#heatscatter(rna_tpm[rna_tpm[,colnum]!=0,colnum]-rna_tpm[rna_tpm[,colnum]!=0,colnum+1], d_raw[rna_tpm[,colnum]!=0,colnum]-d_raw[rna_tpm[,colnum]!=0,colnum+1], log='', main=toString(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson')))
heatscatter(rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1], d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1], log='', main=toString(sum(((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])-(d_raw[used_id_rna_tpm,colnum]-d_raw[used_id_rna_tpm,colnum+1]))^2)/sum((rna_tpm[used_id_rna_tpm,colnum]-rna_tpm[used_id_rna_tpm,colnum+1])^2)) )

abline(0,1,col='red')
print(cor((rna_tpm[,colnum])-rna_tpm[,colnum+1], (d_raw[,colnum])-d_raw[,colnum+1], method='pearson'))
dev.off()





png('check1.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'pknorm', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
plot(rna_tpm[,colnum], d_raw[,colnum], log='')
abline(0,1,col='red')
print(cor((rna_tpm[,colnum]), (d_raw[,colnum])))
dev.off()

png('check2.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'raw', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
plot(rna_tpm[,colnum], d_raw[,colnum], log='')
abline(0,1,col='red')
print(cor((rna_tpm[,colnum]+0.1), (d_raw[,colnum]+0.1)))
dev.off()



png('check1r.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'pknorm', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
plot(rna_tpm[,colnum], d_raw[,colnum+add], log='')
abline(0,1,col='red')
print(cor((rna_tpm[,colnum]+0.1), (d_raw[,colnum]+0.1)))
dev.off()

png('check2r.png')
d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', 'raw', '.txt', sep=''), header=F))+0.1)
d_raw = (d_raw - mean(d_raw))/sd(d_raw)
plot(rna_tpm[,colnum], d_raw[,colnum+add], log='')
abline(0,1,col='red')
print(cor((rna_tpm[,colnum]+0.1), (d_raw[,colnum]+0.1)))
dev.off()

#for (i in c(1:11)){
#colnum = 1
#add=i
#png(paste('check3.', toString(i), '.png', sep=''))
#d_raw = read.table(paste('tss_h3k4me3.pcsorted.', 'raw', '.txt', sep=''), header=F)
#plot(rna_tpm[,colnum], rna_tpm[,colnum+add], log='xy')
#print(cor(log2(rna_tpm[,colnum]), log2(rna_tpm[,colnum])))
#abline(0,1,col='red')
#dev.off()
#}

methods = c('raw', 'trnorm', 'manorm', 'znorm', 'qtnorm', 'pknorm')
k = 10
my.cols = rev(brewer.pal(k, "RdYlBu"))
small_num = 0.01
upper_lim = max((cbind(log2(rna_tpm[,]+small_num), log2(d_raw[,]+small_num))))
lower_lim = min((cbind(log2(rna_tpm[,]+small_num), log2(d_raw[,]+small_num))))

upper_lim = 2
lower_lim = 0


for (i in c(1:11)){
print(paste('sc', toString(i), '.compare.png', sep=''))
png(paste('sc', toString(i), '.compare.png', sep=''), width=1200, height=800)
par(mfrow=c(2,3))
for (m in methods){
d_raw = read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F)
rna_tpm_dif = (log2(rna_tpm[used_id_rna_tpm,i]+small_num))
chip_dif = (log2(d_raw[used_id_rna_tpm,i]+small_num))
dif_cor_1 = cor(rna_tpm_dif, chip_dif, method='pearson')
dif_cor_2 = cor(rna_tpm_dif, chip_dif, method='spearman')
plot(main=paste(m, ': pr: ', toString(round(dif_cor_1, digits=3)), '; sr: ', toString(round(dif_cor_2, digits=3)), sep=''), rna_tpm_dif, chip_dif, log='', xlim=c(lower_lim,upper_lim), ylim=c(lower_lim,upper_lim), pch=16)
z <- kde2d(rna_tpm_dif, chip_dif, n=20)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
abline(0,1,col='red')
abline(h=0, col='blue')
abline(v=0, col='blue')
}
dev.off()
}




###### get rna_dif vs chip-dif 6 plots
methods = c('raw', 'trnorm', 'manorm', 'znorm', 'qtnorm', 'pknorm')
k = 10
cor_matrix = c()
cor_matrix_bg=c()
rna_tpm_max = apply(log2(rna_tpm+small_num), 1, max)
rna_tpm_max_lim = 3
id = c()
y=0
for (i in c(1:11)){
	for (j in c(1:11)){
		all_cor = c()
		all_cor_bg = c()
		x=0
		if (i < j){
			y=y+1
			id[y] = paste(toString(i), '_', toString(j), sep='')

			print(paste('sc', toString(i), '_', toString(j), '.dif.png', sep=''))
			png(paste('sc', toString(i), '_', toString(j), '.dif.png', sep=''), width=1200, height=800)
			par(mfrow=c(2,3))
			for (m in methods){
				d_raw = read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F)
				d_raw_bg = read.table(paste('tss_h3k4me3.pcsorted.', m, '.1000kb.txt', sep=''), header=F)
				#rna_tpm_dif = scale(log2(rna_tpm[rna_tpm_max>rna_tpm_max_lim,i]+small_num)-log2(rna_tpm[rna_tpm_max>rna_tpm_max_lim,j]+small_num), center = FALSE)
				#chip_dif = scale(log2(d_raw[rna_tpm_max>rna_tpm_max_lim,i]+small_num)-log2(d_raw[rna_tpm_max>rna_tpm_max_lim,j]+small_num), center = FALSE)
				#chip_dif_bg = scale(log2(d_raw_bg[rna_tpm_max>rna_tpm_max_lim,i]+small_num)-log2(d_raw_bg[rna_tpm_max>rna_tpm_max_lim,j]+small_num), center = FALSE)
				rna_tpm_dif = log2(rna_tpm[used_id_rna_tpm,i]+small_num)-log2(rna_tpm[used_id_rna_tpm,j]+small_num)
				chip_dif = log2(d_raw[used_id_rna_tpm,i]+small_num)-log2(d_raw[used_id_rna_tpm,j]+small_num)
				chip_dif_bg = log2(d_raw_bg[used_id_rna_tpm,i]+small_num)-log2(d_raw_bg[used_id_rna_tpm,j]+small_num)
				dif_cor_1 = cor(rna_tpm_dif, chip_dif, method='pearson')
				dif_cor_1_bg = cor(rna_tpm_dif, chip_dif_bg, method='pearson')
				dif_cor_2 = cor(rna_tpm_dif, chip_dif, method='spearman')
				plot(main=paste(m, ': pr: ', toString(round(dif_cor_1, digits=3)), '; sr: ', toString(round(dif_cor_2, digits=3)), sep=''), rna_tpm_dif, chip_dif, log='', xlim=c(-15,15), ylim=c(-15,15), pch=16)
				z <- kde2d(rna_tpm_dif, chip_dif, n=20)
				contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
				abline(0,1,col='red')
				abline(h=0, col='blue')
				abline(v=0, col='blue')
				x=x+1
				all_cor[x] = dif_cor_1
				all_cor_bg[x] = dif_cor_1_bg
			}
			dev.off()
			cor_matrix = rbind(cor_matrix, all_cor)
			cor_matrix_bg = rbind(cor_matrix_bg, all_cor_bg)
		}
	}
}

cbind(id,cor_matrix[,6]-apply(cor_matrix, 1, mean))[order(cor_matrix[,6]-apply(cor_matrix, 1, mean)),]

###### get rna_dif vs chip-dif distribution 
for (i in c(1:6)){
	kl_dist_cor_dist_bg = kl.dist(density(cor_matrix[,i], bw=bw_used)$y, density(cor_matrix_bg[,i], bw=bw_used)$y)$D2
	paired_t = t.test(cor_matrix[,i], cor_matrix_bg[,i], paired=TRUE, alternative = 'greater')
	paired_t_statistic = paired_t$statistic
	#-log10(ks.test(cor_dif_matrix_shuffle[,i], cor_dif_matrix[,i], alternative='greater')$p)#
	#print(ks.test(cor_dif_matrix_bg[,i], cor_dif_matrix[,i], alternative='greater'))
	#print(ks.test(cor_dif_matrix_bg[,i], cor_dif_matrix[,i], alternative='greater')$statistic)
	png(paste('tss_h3k4me3.pcsorted.difcor.', methods[i], '.png', sep=''))
	plot(density(cor_matrix[,i], bw=bw_used), col='green', main=paste('paired_t_statistic = ', toString(round(paired_t_statistic, digits=3)), sep=''), ylim=c(0,6))
	#lines(density(cor_dif_matrix_shuffle[,i], bw=bw_used), col='black')
	lines(density(cor_matrix_bg[,i], bw=bw_used), col='blue')
	dev.off()
}

png('test.png')
boxplot((as.matrix(cor_matrix)[,order(colMeans(cor_matrix))]), use.cols = TRUE, xaxt='n')
axis(1, at=1:6, labels=methods[order(colMeans(cor_matrix))])
dev.off()


png('test_5.png')
boxplot((as.matrix(cor_matrix[,-4])[,order(colMeans(cor_matrix[,-4]))]), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=methods[-4][order(colMeans(cor_matrix[,-4]))])
dev.off()

###### get rna_dif vs chip-dif boxplot 6 plots 
###### get table for ggplot
cor_matrix_table = c()
for (i in c(1:dim(cor_matrix)[2])){
	paired_t = t.test(cor_matrix[,i], cor_matrix_bg[,i], paired=TRUE, alternative = 'greater')
	paired_t_statistic = paired_t$statistic
	cor_matrix_table=rbind(cor_matrix_table, cbind(cor_matrix[,i], rep(paste(toString(i), '_', methods[i], ': ', toString(round(paired_t_statistic, digits=3)), sep=''), dim(cor_matrix)[1]), rep('tss',dim(cor_matrix)[1]) ))
}
for (i in c(1:dim(cor_matrix_bg)[2])){
	paired_t = t.test(cor_matrix[,i], cor_matrix_bg[,i], paired=TRUE, alternative = 'greater')
	paired_t_statistic = paired_t$statistic
	cor_matrix_table=rbind(cor_matrix_table, cbind(cor_matrix_bg[,i], rep(paste(toString(i), '_', methods[i], ': ', toString(round(paired_t_statistic, digits=3)), sep=''), dim(cor_matrix)[1]), rep('1000kb',dim(cor_matrix)[1])  ))
}

cor_matrix_table = as.data.frame(cor_matrix_table)
colnames(cor_matrix_table) = c('cor', 'method', 'fg_bg')
cor_matrix_table[,1] = apply(cor_matrix_table, 1, function(x) x[1]=as.numeric(x[1]))

png('test1.png', width=1200, height=800)
par(mfrow=c(1,6))
p = ggplot(data = cor_matrix_table, aes(x=method, y=cor)) 
p = p + geom_boxplot(aes(fill = fg_bg))
p = p + geom_point(aes(y=cor, group=fg_bg), position = position_dodge(width=0.75))
p = p + stat_compare_means(aes(group = fg_bg), label = "p.format", paired = TRUE, method = "t.test")
p = p + facet_wrap( ~ method, scales="free")
p = p + xlab("Methods") + ylab("Pearson Correlation") + ggtitle("Pearson Correlation between different methods")
p = p + guides(fill=guide_legend(title="position"))
p = p + coord_cartesian(ylim = c(min(cor_matrix_table[,1]), max(cor_matrix_table[,1])))
p
dev.off()

library(mixtools)
library(mclust)
mixmdl2 = normalmixEM(log2(rna_tpm)[is.finite(log2(rna_tpm))], k=2)
png('rna_tpm_hist.png')
plot(mixmdl2, which=2)
dev.off()















'''

methods = c('pknorm')
for (m in methods){
	print(m)
	set.seed(1)
	bw_used = 0.1
	cor_method = 'pearson'
	small_num = 0.1
	tpm_lim=2
	shuffle_id = sample(dim(rna_tpm)[1],dim(rna_tpm)[1])
	d_raw = read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F)
	cor_0 = apply(cbind(rna_tpm, d_raw)[rna_tpm_max>=tpm_lim,], 1, function(x) cor(log2(x[1:11]+small_num),log2(x[12:22]+small_num), method=cor_method))
	cor_0_shuffle = apply(cbind(rna_tpm, d_raw[shuffle_id,])[rna_tpm_max>=tpm_lim,], 1, function(x) cor(log2(x[1:11]+small_num),log2(x[12:22]+small_num), method=cor_method))
	d_raw_bg = read.table(paste('tss_h3k4me3.pcsorted.', m, '.1000kb.txt', sep=''), header=F)
	cor_0_bg = apply(cbind(rna_tpm, d_raw_bg)[rna_tpm_max>=tpm_lim,], 1, function(x) cor(log2(x[1:11]+small_num),log2(x[12:22]+small_num), method=cor_method))
	kl_dist_bg = kl.dist(density(cor_0[!is.na(cor_0)], bw=bw_used)$y, density(cor_0_bg[!is.na(cor_0_bg)], bw=bw_used)$y)$D2
	kl_dist_shuffle = kl.dist(density(cor_0[!is.na(cor_0)], bw=bw_used)$y, density(cor_0_shuffle[!is.na(cor_0_shuffle)], bw=bw_used)$y)$D2
}

c_pk = c()
c_qt = c()


for (i in c(1:11)){
	m_b='pknorm_density'
	m_c='qtnorm'
	cor_method = 'spearman'
	d_raw_pk = read.table(paste('tss_h3k4me3.pcsorted.', m_b, '.txt', sep=''), header=F)
	d_raw_qt = read.table(paste('tss_h3k4me3.pcsorted.', m_c, '.txt', sep=''), header=F)
	used_id = (abs(log2(rna_tpm[,i]+0.01)-log2(rna_tpm[,i+1]+0.01)) > 2) 
	print(sum(used_id))
	a=(log2(as.vector((rna_tpm[used_id,i]))+0.01)-log2(as.vector((rna_tpm[used_id,i+1]))+0.01))
	b=(log2(as.vector(as.matrix((d_raw_pk[used_id,i])))+0.01)-log2(as.vector(as.matrix((d_raw_pk[used_id,i+1])))+0.01))
	c=(log2(as.vector(as.matrix((d_raw_qt[used_id,i])))+0.01)-log2(as.vector(as.matrix((d_raw_qt[used_id,i+1])))+0.01))
	par(mfrow=c(1,2))
	plot((a),(b), log='', xlim=c(min(cbind(a,b)), max(cbind(a,b)))+0.01, ylim=c(min(cbind(a,b)), max(cbind(a,b)))+0.01, main=toString(sum(a*b>=0)/sum(a*b!='1')))#main=toString(cor(a,b,method=cor_method)))
	abline(0,1,col='red')
	plot((a),(c), log='', xlim=c(min(cbind(a,c)), max(cbind(a,c)))+0.01, ylim=c(min(cbind(a,c)), max(cbind(a,c)))+0.01, main=toString(sum(a*c>=0)/sum(a*c!='1')))#main=toString(cor(a,c,method=cor_method)))
	abline(0,1,col='red')
	print(i)
	print(cor(a,b))
	print(cor(a,c))
	c_pk[i] = cor(a,b, method=cor_method)
	c_qt[i] = cor(a,c, method=cor_method)
}


cbind(c_pk, c_qt)

plot(rna_tpm[rna_tpm_max>=tpm_lim,3]/rna_tpm[rna_tpm_max>=tpm_lim,2], d_raw[rna_tpm_max>=tpm_lim,3]/d_raw[rna_tpm_max>=tpm_lim,2], log='xy', xlim=c(0.0001,10000), ylim=c(0.0001,10000),pch=16)
abline(0,1, col='red')

methods = c('raw', 'pknorm', 'loessnorm', 'qtnorm', 'trnorm', 'manorm')
library(MASS)
library(RColorBrewer)
k <- 10
my.cols <- rev(brewer.pal(k, "RdYlBu"))

cor_dif_matrix = c()
cor_dif_matrix_bg = c()
cor_dif_matrix_shuffle = c()
tpm_lim=5
used_gene_number = c()
count_id = 0
cor_method = 'pearson'
for (m in methods){
	print(m)
	d_raw = read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F)
	d_raw_bg = read.table(paste('tss_h3k4me3.pcsorted.', m, '.1000kb.txt', sep=''), header=F)
	cor_dif_col = c()
	cor_dif_col_bg = c()
	cor_dif_col_shuffle = c()
	add_small_num = 0.01
	for (i in c(1:10)){
		for (j in c((i+1):11)){
			print(paste(toString(i), '_', toString(j), sep=''))
			png(paste('tss_h3k4me3.pcsorted.dif.', m, toString(i), '_', toString(j), '.png', sep=''))
			rna_tpm_pair = apply(rna_tpm[,c(i,j)],1, max)
			rna_log1 = log2(rna_tpm[rna_tpm_pair>=tpm_lim,i]+add_small_num)
			rna_log2 = log2(rna_tpm[rna_tpm_pair>=tpm_lim,j]+add_small_num)
			chip_log1 = log2(d_raw[rna_tpm_pair>=tpm_lim,i]+add_small_num)
			chip_log2 = log2(d_raw[rna_tpm_pair>=tpm_lim,j]+add_small_num)
			chip_log1_bg = log2(d_raw_bg[rna_tpm_pair>=tpm_lim,i]+add_small_num)
			chip_log2_bg = log2(d_raw_bg[rna_tpm_pair>=tpm_lim,j]+add_small_num)
			rna_dif = scale(rna_log1-rna_log2, center = FALSE)
			chip_dif = scale(chip_log1-chip_log2, center = FALSE)
			chip_dif_bg = scale(chip_log1_bg-chip_log2_bg, center = FALSE)
			cor_dif = cor(rna_dif, chip_dif, method=cor_method)
			count_id = count_id+1
			used_gene_number[count_id] = sum(rna_tpm_pair>=tpm_lim)
			set.seed(2018)
			shuffle_id_dif = sample(length(rna_dif), length(rna_dif))
			cor_dif_shuffle = cor(rna_dif, chip_dif[shuffle_id_dif], method=cor_method)
			cor_dif_bg = cor(rna_dif, chip_dif_bg, method=cor_method)
			cor_dif_col = rbind(cor_dif_col, cor_dif)
			cor_dif_col_bg = rbind(cor_dif_col_bg, cor_dif_bg)
			cor_dif_col_shuffle = rbind(cor_dif_col_shuffle, cor_dif_shuffle)
			z <- kde2d(rna_dif, chip_dif, n=50)
			plot(scale(rna_log1-rna_log2, center = FALSE), scale(chip_log1-chip_log2, center = FALSE), log='', xlim=c(-5,5), ylim=c(-5,5),pch=16)
			contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
			abline(0,1, col='red')
			abline(h=0, col='blue')
			abline(v=0, col='blue')
			dev.off()
		}
	}
	cor_dif_matrix = cbind(cor_dif_matrix, cor_dif_col)
	cor_dif_matrix_bg = cbind(cor_dif_matrix_bg, cor_dif_col_bg)
	cor_dif_matrix_shuffle = cbind(cor_dif_matrix_shuffle, cor_dif_col_shuffle)
}

library(mixtools)
library(mclust)
mixmdl2 = normalmixEM(log2(rna_tpm)[is.finite(log2(rna_tpm))], k=2)
png('rna_tpm_hist.png')
plot(mixmdl2, which=2)
dev.off()

png('used_gene_number_hist.png')
hist((used_gene_number), breaks=20)
dev.off()

methods = c('raw', 'pknorm', 'loessnorm', 'qtnorm', 'trnorm', 'manorm')

for (i in c(1:6)){
	kl_dist_cor_dist_bg = kl.dist(density(cor_dif_matrix[,i], bw=bw_used)$y, density(cor_dif_matrix_bg[,i], bw=bw_used)$y)$D2
	kl_dist_cor_dist_shuffle = kl.dist(density(cor_dif_matrix[,i], bw=bw_used)$y, density(cor_dif_matrix_shuffle[,i], bw=bw_used)$y)$D2
	paired_t = t.test(cor_dif_matrix[,i], cor_dif_matrix_bg[,i], paired=TRUE, alternative = 'greater')
	paired_t_statistic = paired_t$statistic
	#-log10(ks.test(cor_dif_matrix_shuffle[,i], cor_dif_matrix[,i], alternative='greater')$p)#
	#print(ks.test(cor_dif_matrix_bg[,i], cor_dif_matrix[,i], alternative='greater'))
	#print(ks.test(cor_dif_matrix_bg[,i], cor_dif_matrix[,i], alternative='greater')$statistic)
	png(paste('tss_h3k4me3.pcsorted.difcor.', methods[i], '.png', sep=''))
	plot(density(cor_dif_matrix[,i], bw=bw_used), col='green', main=paste('paired_t_statistic = ', toString(round(paired_t_statistic, digits=3)), sep=''), ylim=c(0,6))
	#lines(density(cor_dif_matrix_shuffle[,i], bw=bw_used), col='black')
	lines(density(cor_dif_matrix_bg[,i], bw=bw_used), col='blue')
	dev.off()
}



library(ggpubr)

png('test.png')
FRiP_dif_flatt = as.vector(cor_dif_matrix)
method_name = rep(methods, dim(cor_dif_matrix)[2])
FRiP_dif_flatt_matrix = cbind(as.data.frame(FRiP_dif_flatt), as.data.frame(method_name))
FRiP_dif_flatt_matrix = FRiP_dif_flatt_matrix[!is.na(FRiP_dif_flatt_matrix[,1]),]
colnames(FRiP_dif_flatt_matrix) = c('cor_dif' ,'method')
my_comparisons = list( c("PKnorm", "raw"), c("PKnorm", "TRnorm"), c("PKnorm", "MAnorm") )
ggboxplot(FRiP_dif_flatt_matrix, x = "method", y = "cor_dif")+ 
stat_compare_means(comparisons = my_comparisons, paired = TRUE, method='t.test', method.args = list(alternative ="less"), na.rm = TRUE)
dev.off()

png('compare_pk_loess.png')
low_lim = min(cbind(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,3]-cor_dif_matrix_bg[,3]))
up_lim = max(cbind(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,3]-cor_dif_matrix_bg[,3]))
plot(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,3]-cor_dif_matrix_bg[,3], xlim=c(low_lim,up_lim), ylim=c(low_lim,up_lim))
abline(0,1, col='red')
dev.off()


png('compare_pk_qt.png')
low_lim = min(cbind(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,4]-cor_dif_matrix_bg[,4]))
up_lim = max(cbind(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,4]-cor_dif_matrix_bg[,4]))
plot(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,4]-cor_dif_matrix_bg[,4], xlim=c(low_lim,up_lim), ylim=c(low_lim,up_lim))
abline(0,1, col='red')
dev.off()

png('test.png')
boxplot((cor_dif_matrix[,order(colMeans(cor_dif_matrix))]), use.cols = TRUE, xaxt='n')
axis(1, at=1:6, labels=methods[order(colMeans(cor_dif_matrix))])
dev.off()


id=c()
k=0
for (i in c(1:10)){
	for (j in c((i+1):11)){
		k=k+1
		id[k] = paste(toString(i), '_', toString(j), sep='')
	}
}

cbind(id,cor_dif_matrix[,2]-apply(cor_dif_matrix, 1, max))[order(cor_dif_matrix[,2]-apply(cor_dif_matrix, 1, max)),]



apply(cor_dif_matrix, 1, max)


'''
