library(seewave)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

files_list = list.files(pattern = "\\.all5info.txt$")


info_1 = read.table(files_list[1], header=TRUE, sep='\t')
info_type = colnames(info_1)
info_methods = rownames(info_1)


for (i in c(1:length(info_type))){
	info_matrix = c()
	for (j in c(1:length(files_list))){
		info_tmp = read.table(files_list[j], header=TRUE, sep='\t')[i]
		ct = unlist(strsplit(files_list[j], split='[.]'))[1]
		if (info_type[i]!='snr'){ref_sig = info_tmp[1,]}
		if (info_type[i]=='snr'){ref_sig = log2(info_tmp[1,])}
		#for (k in c(2:length(info_methods))){
		for (k in c(2,3,4,5,6)){
			if (info_type[i]!='snr'){
				info_matrix = rbind(info_matrix, c(info_methods[k], (info_tmp[k,])))
			} else{
				info_matrix = rbind(info_matrix, c(info_methods[k], log2(info_tmp[k,])))
			}
		}	
	}
	info_matrix = as.data.frame(info_matrix)
	colnames(info_matrix) = c('method', 'sig')
	info_matrix$method = factor(info_matrix$method, levels = c('sig2_r','sig2_tr','sig2_ma','sig2_qt','sig2_pk'),ordered = TRUE)
	#print(info_matrix)
	info_matrix[,2] = apply(info_matrix, 1, function(x) as.numeric(x[2]))
	pdf(paste(info_type[i], '.box.pdf', sep=''))#, width=500, height=500)
	p = ggplot(data = info_matrix, aes(x=method, y=sig)) 
	p = p + geom_boxplot(aes(fill = method))
	p = p + geom_point(aes(y=sig, group=method), position = position_dodge(width=0.75))
	p = p + geom_hline(yintercept = ref_sig, color="black", linetype="dashed")
	p = p + scale_fill_manual(values=c("white", "white", "white",'white','white')) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
	#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
	plot(p)
	dev.off()
}










