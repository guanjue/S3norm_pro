total_bin_num = 200000
pk_bins_num = 20000
dif_pk_bins_num = 4000
dif_fc = 4

all_bins = rep(0, total_bin_num)
### rep1 spe pk
all_bins[1:dif_pk_bins_num]=1
### rep2 spe pk
all_bins[(dif_pk_bins_num+1):(dif_pk_bins_num*2)]=2
### common pk
all_bins[(dif_pk_bins_num*2+1):pk_bins_num]=3


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
p = ggplot(data = cor_matrix_table, aes(x=method, y=cor)) 
p = p + geom_boxplot(aes(fill = fg_bg))
p = p + geom_point(aes(y=cor, group=fg_bg), position = position_dodge(width=0.75))
p = p + facet_wrap( ~ method, scales="free")
p = p + xlab("Methods") + ylab("Pearson Correlation") + ggtitle("Pearson Correlation between different methods")
p = p + guides(fill=guide_legend(title="position"))
p = p + coord_cartesian(ylim = c(min(cor_matrix_table[,1]), max(cor_matrix_table[,1])))
p
dev.off()


