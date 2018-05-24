library(ggpubr)

files = list.files(pattern = "\\.txt$")

jaccard_index = c()
FRiP_dif = c()
adjusted_rand_index = c()
FRiP_JI_all = c()
FRuP_JI_all = c()
FRiP_FRuP = c()
ct_id = c()
FRiP = c()
i=1
for (file in files){
	print(file)
	file_name_split = unlist(strsplit(file, "_"))	
	info = read.table(file, header=T)
	info = info[c(1,4,5,2,3),]
	print(info)
	if (dim(info)[2]!=2000){
		FRiP_dif = cbind(FRiP_dif, abs(info[,2]-info[,1]))### FRiP_dif
		FRiP = cbind(FRiP, info[,1])
		jaccard_index = cbind(jaccard_index, (info[,5]))### jaccard_index
		adjusted_rand_index = cbind(adjusted_rand_index, (info[,6]))### adjusted_rand_index
		FRiP_JI_all = cbind(FRiP_JI_all, (info[,15]))### FRiP_JI_all
		FRuP_JI_all = cbind(FRuP_JI_all, (info[,18]))### FRuP_JI_all
		FRiP_FRuP = cbind(FRiP_FRuP, (info[,15]/info[,18]))### FRiP_FRuP
		ct_id[i] = paste(file_name_split[1], file_name_split[2], file_name_split[3], sep='_')
		i = i+1
	}
	
}

rownames(jaccard_index) = rownames(info)
colnames(jaccard_index) = ct_id
jaccard_index

png('test.png')
boxplot((cor_dif_matrix), use.cols = TRUE)
dev.off()


axis(1, at=1:5, labels=rownames(jaccard_index))

pdf('FRiP_dif.pdf')
boxplot(t(FRiP_dif), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRiP_dif_sig.pdf')
#boxplot(t(FRiP_dif), use.cols = TRUE, xaxt='n')
#axis(1, at=1:5, labels=rownames(jaccard_index))
FRiP_dif_flatt = as.vector(FRiP_dif)
method_name = rep(rownames(jaccard_index), dim(FRiP_dif)[2])
FRiP_dif_flatt_matrix = cbind(as.data.frame(FRiP_dif_flatt), as.data.frame(method_name))
colnames(FRiP_dif_flatt_matrix) = c('FRiP_dif' ,'method')
my_comparisons = list( c("PKnorm", "raw"), c("PKnorm", "TRnorm"), c("PKnorm", "MAnorm"), c("PKnorm", "QTnorm") )
ggboxplot(FRiP_dif_flatt_matrix, x = "method", y = "FRiP_dif")+ 
stat_compare_means(comparisons = my_comparisons, paired = TRUE, method='t.test', method.args = list(alternative ="less"))
dev.off()


pdf('FRiP_dif_sig_4.pdf')
#boxplot(t(FRiP_dif), use.cols = TRUE, xaxt='n')
#axis(1, at=1:5, labels=rownames(jaccard_index))
FRiP_dif_flatt = as.vector(FRiP_dif)
method_name = rep(rownames(jaccard_index), dim(FRiP_dif)[2])
FRiP_dif_flatt_matrix = cbind(as.data.frame(FRiP_dif_flatt), as.data.frame(method_name))
FRiP_dif_flatt_matrix = FRiP_dif_flatt_matrix[FRiP_dif_flatt_matrix[,2]!='QTnorm',]
colnames(FRiP_dif_flatt_matrix) = c('FRiP_dif' ,'method')
my_comparisons = list( c("PKnorm", "raw"), c("PKnorm", "TRnorm"), c("PKnorm", "MAnorm") )
ggboxplot(FRiP_dif_flatt_matrix, x = "method", y = "FRiP_dif")+ 
stat_compare_means(comparisons = my_comparisons, paired = TRUE, method='t.test', method.args = list(alternative ="less"))
dev.off()


pdf('FRiP_sig.pdf')
#boxplot(t(FRiP_dif), use.cols = TRUE, xaxt='n')
#axis(1, at=1:5, labels=rownames(jaccard_index))
FRiP_flatt = as.vector(FRiP)
method_name = rep(rownames(jaccard_index), dim(FRiP)[2])
FRiP_flatt_matrix = cbind(as.data.frame(FRiP_flatt), as.data.frame(method_name))
colnames(FRiP_flatt_matrix) = c('FRiP' ,'method')
my_comparisons = list( c("PKnorm", "raw"), c("PKnorm", "TRnorm"), c("PKnorm", "MAnorm"), c("PKnorm", "QTnorm") )
ggboxplot(FRiP_flatt_matrix, x = "method", y = "FRiP")+ 
stat_compare_means(comparisons = my_comparisons, paired = TRUE, method='t.test', method.args = list(alternative ="greater"))
dev.off()


pdf('jaccard_index.pdf')
boxplot(t(jaccard_index), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('adjusted_rand_index.pdf')
boxplot(t(adjusted_rand_index), use.cols = TRUE, xaxt='n', ylim=c(0,1))
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRiP_JI_all.pdf')
boxplot(t(FRiP_JI_all), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRuP_JI_all.pdf')
boxplot(t(FRuP_JI_all), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRiP_FRuP.pdf')
boxplot(t(FRiP_FRuP), use.cols = TRUE, xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()



pdf('FRiP_dif_lines.pdf')
plot((FRiP_dif[,1]), type='l', xaxt='n', ylim=c(min(FRiP_dif), max(FRiP_dif)))
for (i in seq(2,dim(FRiP_dif)[2])){
	lines((FRiP_dif[,i]))
}
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()


pdf('adjusted_rand_index_lines.pdf')
plot((adjusted_rand_index[,1]), type='l', xaxt='n', ylim=c(min(adjusted_rand_index), max(adjusted_rand_index)))
for (i in seq(2,dim(adjusted_rand_index)[2])){
	lines((adjusted_rand_index[,i]))
}
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()





pdf('jaccard_index_ER4_98.pdf')
plot((jaccard_index[,3]), type='p', xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRiP_JI_all_ER4_98.pdf')
plot((FRiP_JI_all[,3]), type='p', xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRuP_JI_all_ER4_98.pdf')
plot((FRuP_JI_all[,3]), type='p', xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()

pdf('FRiP_FRuP_ER4_98.pdf')
plot((FRiP_FRuP[,3]), type='p', xaxt='n')
axis(1, at=1:5, labels=rownames(jaccard_index))
dev.off()









boxplot(t(jaccard_index), use.cols = TRUE)

jaccard_index_rank = apply(jaccard_index,2,function(x) rank(x))

jaccard_index_dif = apply(jaccard_index,1,function(x) x-jaccard_index[1,])

boxplot(jaccard_index_dif, use.cols = TRUE)

lines_plot = t(jaccard_index_dif)

lines_plot = (jaccard_index)
plot(lines_plot[,1], type='l', ylim=c(min(lines_plot), max(lines_plot)), xaxt='n')

for (i in c(2:dim(lines_plot)[2])){
	lines(lines_plot[,i])
}

axis(1, at=1:5, labels=rownames(jaccard_index))


boxplot(t(jaccard_index_rank), use.cols = TRUE)

boxplot(t(jaccard_index), use.cols = TRUE)


