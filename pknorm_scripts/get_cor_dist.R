library(seewave)

data0 =read.table('rna_rpk.pcsorted.txt', header=F)

total_rpk = colSums(data0[,-1])
rna_tpm = t(apply(data0[,-1], 1, function(x) x/total_rpk*1000000))
rna_tpm_max = apply(rna_tpm, 1, max)
rna_tpm_min = apply(rna_tpm, 1, min)

methods = c('raw', 'pknorm', 'qtnorm', 'trnorm', 'manorm')

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
	png(paste('tss_h3k4me3.pcsorted.', m, '.png', sep=''))
	plot(density(cor_0[!is.na(cor_0)], bw=bw_used), col='green', main=paste('KL-dist = ', toString(kl_dist_bg), ' ', toString(kl_dist_shuffle), sep=''), ylim=c(0,1.2))
	lines(density(cor_0_bg[!is.na(cor_0_bg)], bw=bw_used), col='black')
	lines(density(cor_0_shuffle[!is.na(cor_0_shuffle)], bw=bw_used), col='blue')
	dev.off()
}

methods = c('trnorm')

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

methods = c('raw', 'pknorm', 'qtnorm', 'trnorm', 'manorm')
library(MASS)
library(RColorBrewer)
k <- 10
my.cols <- rev(brewer.pal(k, "RdYlBu"))

cor_dif_matrix = c()
cor_dif_matrix_bg = c()
cor_dif_matrix_shuffle = c()
tpm_lim=1
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
			rna_dif = (rna_log1-rna_log2)
			chip_dif = (chip_log1-chip_log2)
			chip_dif_bg = (chip_log1_bg-chip_log2_bg)
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
			plot(rna_log1-rna_log2, chip_log1-chip_log2, log='', xlim=c(-10,10), ylim=c(-10,10),pch=16)
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


for (i in c(1:5)){
	kl_dist_cor_dist_bg = kl.dist(density(cor_dif_matrix[,i], bw=bw_used)$y, density(cor_dif_matrix_bg[,i], bw=bw_used)$y)$D2
	kl_dist_cor_dist_shuffle = kl.dist(density(cor_dif_matrix[,i], bw=bw_used)$y, density(cor_dif_matrix_shuffle[,i], bw=bw_used)$y)$D2
	#-log10(ks.test(cor_dif_matrix_shuffle[,i], cor_dif_matrix[,i], alternative='greater')$p)#
	#print(ks.test(cor_dif_matrix_bg[,i], cor_dif_matrix[,i], alternative='greater'))
	#print(ks.test(cor_dif_matrix_bg[,i], cor_dif_matrix[,i], alternative='greater')$statistic)
	png(paste('tss_h3k4me3.pcsorted.difcor.', methods[i], '.png', sep=''))
	plot(density(cor_dif_matrix[,i], bw=bw_used), col='green', main=paste('KL-dist = ', toString(round(kl_dist_cor_dist_bg, digits=3)), sep=''), ylim=c(0,6))
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

png('compare.png')
low_lim = min(cbind(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,3]-cor_dif_matrix_bg[,3]))
up_lim = max(cbind(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,3]-cor_dif_matrix_bg[,3]))
plot(cor_dif_matrix[,2]-cor_dif_matrix_bg[,2], cor_dif_matrix[,3]-cor_dif_matrix_bg[,3], xlim=c(low_lim,up_lim), ylim=c(low_lim,up_lim))
abline(0,1, col='red')
dev.off()

png('test.png')
boxplot((cor_dif_matrix), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=methods)
dev.off()




