
data0 =read.table('rna_rpk.pcsorted.txt', header=F)
#data0 =read.table('rna_rpk.pcsorted.all12.txt', header=F)

total_rpk = colSums(data0[,-1])

rna_tpm = t(apply(data0[,-1], 1, function(x) x/total_rpk*10000000))

rna_tpm_max = apply(rna_tpm, 1, max)
tpm_lim=5
pdf('hist_rna_tpm_max.pdf')
hist(log2(rna_tpm_max), breaks=50)
abline(v=tpm_lim, col='red', lwd=1.5, lty=3)
dev.off()

used_id_rna_tpm = ( (log2(rna_tpm_max+0.1)>tpm_lim)  ) >0
print(sum(used_id_rna_tpm))


rna_tpm = as.matrix(log10(rna_tpm+0.1))
rna_tpm = (rna_tpm-mean(rna_tpm))/sd(rna_tpm)

library(LSD)
ct_list = c('CFU_E_ad', 'CMP', 'ERY_ad', 'GMP', 'MK_imm_ad', 'LSK_BM', 'MEP', 'MONO_BM', 'NEU', 'ER4', 'G1E')
methods_check = c('rcnorm', 'poisnorm', 'rcznorm', 'raw', 'trnorm', 'qtnorm', 'manorm', 'pknorm')
methods_check_name = c('RC', 'POIS', 'Z', 'NB', 'NB-TRnorm', 'NB-QTnorm', 'NB-MAnorm', 'NB-PKnorm')

methods_check_shuf = c('rcnorm_shuf', 'poisnorm_shuf', 'rcznorm_shuf', 'raw_shuf', 'trnorm_shuf', 'qtnorm_shuf', 'manorm_shuf', 'pknorm_shuf')
methods_mix = c('rcnorm', 'rcnorm_shuf', 'poisnorm', 'poisnorm_shuf', 'rcznorm', 'rcznorm_shuf', 'raw', 'raw_shuf', 'trnorm', 'trnorm_shuf', 'qtnorm', 'qtnorm_shuf', 'manorm', 'manorm_shuf', 'pknorm', 'pknorm_shuf')

r_mat = c()
r_mat_sp = c()
r_mat_shuf = c()
p_mat = c()
p_mat_shuf = c()
ct_pair_vec = c()

set.seed(2018)
shuffle_id1 = sample(1:11,11,replace=F)
shuffle_id2 = sample(1:11,11,replace=F)

for (i in c(1:11)){
	for (j in c(1:11)){
		if (i<j){
			k=0
			ct_pair_vec = rbind(ct_pair_vec, paste(ct_list[i], '-', ct_list[j], sep=''))
			r_vec = c()
			r_vec_sp = c()
			r_vec_shuf = c()
			p_vec = c()
			p_vec_shuf = c()

			m = 'rcnorm'
			d_raw0 = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))+0.1)
			d_raw0 = (d_raw0 - mean(d_raw0))/sd(d_raw0)
			rna_dif0 = rna_tpm[used_id_rna_tpm,i]-rna_tpm[used_id_rna_tpm,j]
			hist_dif0 = d_raw0[used_id_rna_tpm,i]-d_raw0[used_id_rna_tpm,j]

			rss0 = (sum((rna_dif0 - hist_dif0)^2))
			for (m in methods_check){
				k = k+1
				png(paste('tss_h3k4me3.pcsorted.check.', m, '.', ct_list[i], '_', ct_list[j], '.png', sep=''))
				if (m != 'rcznorm'){
					d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))+0.1)
				} else {
					d_raw = as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))
				}
				d_raw = (d_raw - mean(d_raw))/sd(d_raw)
				rna_dif = rna_tpm[used_id_rna_tpm,i]-rna_tpm[used_id_rna_tpm,j]
				hist_dif = d_raw[used_id_rna_tpm,i]-d_raw[used_id_rna_tpm,j]
				hist_dif_shuffle = d_raw[used_id_rna_tpm,shuffle_id1[i]]-d_raw[used_id_rna_tpm,shuffle_id1[j]]
				perform = (sum((rna_dif - hist_dif)^2)) / rss0
				perform_shuf = (sum((rna_dif - hist_dif_shuffle)^2)) / rss0
				cor0 = cor(rna_dif, hist_dif, method = 'pearson')
				cor0_sp = cor(rna_dif, hist_dif, method = 'spearman')
				cor0_shuf = cor(rna_dif, hist_dif_shuffle, method = 'pearson')
				r_vec[k] = cor0
				r_vec_sp[k] = cor0_sp
				r_vec_shuf[k] = cor0_shuf
				p_vec[k] = perform
				p_vec_shuf[k] = perform_shuf
				heatscatter(rna_tpm[used_id_rna_tpm,i]-rna_tpm[used_id_rna_tpm,j], d_raw[used_id_rna_tpm,i]-d_raw[used_id_rna_tpm,j], log='', main=paste('R = ', toString(round(cor0, digits=3)), '; Ratio of RSS = ', toString(round(perform, digits=3)), sep=''), xlim=c(-3,3), ylim=c(-3,3) )
				abline(0,1,col='red')
				abline(v=0,col='blue')
				abline(h=0,col='blue')
				print(paste(ct_list[i], '_', ct_list[j], ': ', m, ': R = ', toString(round(cor0, digits=3))) )
				print(paste(ct_list[i], '_', ct_list[j], ': ', m, ': RSS = ', toString(round(perform, digits=3))) )
				dev.off()
			}
			r_mat = rbind(r_mat, r_vec)
			r_mat_sp = rbind(r_mat_sp, r_vec_sp)
			r_mat_shuf = rbind(r_mat_shuf, r_vec_shuf)
			p_mat = rbind(p_mat, p_vec)
			p_mat_shuf = rbind(p_mat_shuf, p_vec_shuf)
		}
	}
}



colnames(r_mat) = methods_check
colnames(r_mat_sp) = methods_check
colnames(r_mat_shuf) = methods_check
colnames(p_mat) = methods_check
colnames(p_mat_shuf) = methods_check

rownames(r_mat) = ct_pair_vec
rownames(r_mat_sp) = ct_pair_vec
rownames(r_mat_shuf) = ct_pair_vec
rownames(p_mat) = ct_pair_vec
rownames(p_mat_shuf) = ct_pair_vec

write.table(r_mat, 'r_mat.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(r_mat_sp, 'r_mat_sp.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )

write.table(r_mat_shuf, 'r_mat_shuf.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(p_mat, 'p_mat.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(p_mat_shuf, 'p_mat_shuf.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )

r_mat = read.table('r_mat.txt', header = TRUE)
r_mat_sp = read.table('r_mat_sp.txt', header = TRUE)

r_mat_shuf = read.table('r_mat_shuf.txt', header = TRUE)
p_mat = read.table('p_mat.txt', header = TRUE)
p_mat_shuf = read.table('p_mat_shuf.txt', header = TRUE)


library(seewave)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

cor_mat_for_boxplot = c()
cor_mat_for_boxplot_sp = c()
cor_mat_for_boxplot_shuf = c()
p_mat_for_boxplot = c()
p_mat_for_boxplot_shuf = c()

for (i in c(1:length(methods_check))){
	### r_mat
	r_mat_tmp = r_mat[,i]
	method_tmp = rep(methods_check_name[i], length(r_mat[,i]))
	cor_mat_for_boxplot = rbind(cor_mat_for_boxplot, cbind(method_tmp, r_mat_tmp))
	### r_mat
	r_mat_tmp = r_mat_sp[,i]
	method_tmp = rep(methods_check_name[i], length(r_mat_sp[,i]))
	cor_mat_for_boxplot_sp = rbind(cor_mat_for_boxplot_sp, cbind(method_tmp, r_mat_tmp))
	### r_mat_shuf
	r_mat_tmp = r_mat_shuf[,i]
	method_tmp = rep(paste(methods_check_name[i], '_shuf', sep=''), length(r_mat[,i]))
	cor_mat_for_boxplot_shuf = rbind(cor_mat_for_boxplot_shuf, cbind(method_tmp, r_mat_tmp))
	### p_mat
	p_mat_tmp = p_mat[,i]
	method_tmp = rep(methods_check_name[i], length(p_mat[,i]))
	p_mat_for_boxplot = rbind(p_mat_for_boxplot, cbind(method_tmp, p_mat_tmp))
	### p_mat_shuf
	p_mat_tmp_shuf = p_mat_shuf[,i]
	method_tmp = rep(paste(methods_check_name[i], '_shuf', sep=''), length(p_mat[,i]))
	p_mat_for_boxplot_shuf = rbind(p_mat_for_boxplot_shuf, cbind(method_tmp, p_mat_tmp_shuf))
}

cor_mat_for_boxplot = as.data.frame(cor_mat_for_boxplot)
colnames(cor_mat_for_boxplot) = c('method', 'R')
cor_mat_for_boxplot$method = factor(cor_mat_for_boxplot$method, levels = methods_check_name,ordered = TRUE)
cor_mat_for_boxplot[,2] = apply(cor_mat_for_boxplot, 1, function(x) as.numeric(x[2]))
pdf(paste('R', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = cor_mat_for_boxplot, aes(x=method, y=R)) 
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=R, group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_name))) 
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

cor_mat_for_boxplot_sp = as.data.frame(cor_mat_for_boxplot_sp)
colnames(cor_mat_for_boxplot_sp) = c('method', 'R')
cor_mat_for_boxplot_sp$method = factor(cor_mat_for_boxplot_sp$method, levels = methods_check_name,ordered = TRUE)
cor_mat_for_boxplot_sp[,2] = apply(cor_mat_for_boxplot_sp, 1, function(x) as.numeric(x[2]))
pdf(paste('R_sp', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = cor_mat_for_boxplot_sp, aes(x=method, y=R)) 
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=R, group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_name))) 
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

cor_mat_for_boxplot_shuf = as.data.frame(cor_mat_for_boxplot_shuf)
colnames(cor_mat_for_boxplot_shuf) = c('method', 'R')
cor_mat_for_boxplot_shuf$method = factor(cor_mat_for_boxplot_shuf$method, levels = methods_check_shuf,ordered = TRUE)
cor_mat_for_boxplot_shuf[,2] = apply(cor_mat_for_boxplot_shuf, 1, function(x) as.numeric(x[2]))
pdf(paste('R.shuf', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = cor_mat_for_boxplot_shuf, aes(x=method, y=R)) 
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=R, group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check)))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

p_mat_for_boxplot = as.data.frame(p_mat_for_boxplot)
colnames(p_mat_for_boxplot) = c('method', 'Performance')
p_mat_for_boxplot$method = factor(p_mat_for_boxplot$method, levels = methods_check_name,ordered = TRUE)
p_mat_for_boxplot[,2] = apply(p_mat_for_boxplot, 1, function(x) as.numeric(x[2]))
pdf(paste('performance', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot, aes(x=method, y=(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()


p_mat_for_boxplot_shuf = as.data.frame(p_mat_for_boxplot_shuf)
colnames(p_mat_for_boxplot_shuf) = c('method', 'Performance')
p_mat_for_boxplot_shuf$method = factor(p_mat_for_boxplot_shuf$method, levels = methods_check_shuf,ordered = TRUE)
p_mat_for_boxplot_shuf[,2] = apply(p_mat_for_boxplot_shuf, 1, function(x) as.numeric(x[2]))
pdf(paste('performance.shuf', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot_shuf, aes(x=method, y=(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_shuf))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()


p_mat_for_boxplot_all = as.data.frame(rbind(p_mat_for_boxplot, p_mat_for_boxplot_shuf))
colnames(p_mat_for_boxplot_all) = c('method', 'Performance')
p_mat_for_boxplot_all$method = factor(p_mat_for_boxplot_all$method, levels = methods_mix,ordered = TRUE)
p_mat_for_boxplot_all[,2] = apply(p_mat_for_boxplot_all, 1, function(x) as.numeric(x[2]))
pdf(paste('performance.all', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot_all, aes(x=method, y=1/(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=1/(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
p = p + scale_fill_manual(values=rep(c("deepskyblue", 'gray'), length(methods_mix))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()


