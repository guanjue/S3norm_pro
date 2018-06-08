library(seewave)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

files_list = list.files(pattern = "\\.all5info.txt$")

print(files_list)

info_1 = read.table(files_list[1], header=TRUE, sep='\t')
info_type = colnames(info_1)
info_methods = rownames(info_1)


for (i in c(1:length(info_type))){
	print(i)
	info_matrix = c()
	for (j in c(1:length(files_list))){
		info_tmp = read.table(files_list[j], header=TRUE, sep='\t')[i]
		ct = unlist(strsplit(files_list[j], split='[.]'))[1]
		for (k in c(1:length(info_methods))){
			info_matrix = rbind(info_matrix, c(info_methods[k], info_tmp[k,]))
		}	
	}
	ref_info = info_tmp[1,]
	info_matrix = as.data.frame(info_matrix)
	colnames(info_matrix) = c('method', 'sig')
	info_matrix$method <- factor(info_matrix$method, levels = info_matrix$method)
	info_matrix[,2] = apply(info_matrix, 1, function(x) as.numeric(x[2]))
	if (info_type[i] == 'snr1' | info_type[i] == 'snr2'){
		info_matrix[,2] = log2(info_matrix[,2])
		ref_info = log2(ref_info)
	} 
	print(paste(info_type[i], '.box.pdf', sep=''))
	pdf(paste(info_type[i], '.box.pdf', sep=''))
	p = ggplot(data = info_matrix, aes(x=method, y=sig)) 
	p = p + geom_boxplot(aes(fill = method))
	p = p + geom_point(aes(y=sig, group=method), position = position_dodge(width=0.75))
	#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
	#p = p + geom_hline(yintercept = ref_info, lty=2, col='blue')
	plot(p)
	dev.off()
	if (i <=3){
		info_matrix = c()
		for (j in c(1:length(files_list))){
			info1_tmp = read.table(files_list[j], header=TRUE, sep='\t')[i]
			info2_tmp = read.table(files_list[j], header=TRUE, sep='\t')[i+3]
			ct = unlist(strsplit(files_list[j], split='[.]'))[1]
			for (k in c(1:length(info_methods))){
				info_matrix = rbind(info_matrix, c(info_methods[k], info1_tmp[k,]-info2_tmp[k,]))
			}	
		}
		ref_info = info_tmp[1,]
		info_matrix = as.data.frame(info_matrix)
		colnames(info_matrix) = c('method', 'sig')
		info_matrix$method <- factor(info_matrix$method, levels = info_matrix$method)
		info_matrix[,2] = apply(info_matrix, 1, function(x) as.numeric(x[2]))
		if (info_type[i] == 'snr1'){
			info_matrix[,2] = log2(info_matrix[,2])
			ref_info = log2(ref_info)
		} else if (info_type[i] == 'pk_num1'){
			info_matrix[,2] = log2(info_matrix[,2])
			ref_info = log2(ref_info)
		}
		print(paste(info_type[i], '.dif.box.pdf', sep=''))
		pdf(paste(info_type[i], '.dif.box.pdf', sep=''))
		p = ggplot(data = info_matrix, aes(x=method, y=sig)) 
		p = p + geom_boxplot(aes(fill = method))
		p = p + geom_point(aes(y=sig, group=method), position = position_dodge(width=0.75))
		#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
		#p = p + geom_hline(yintercept = ref_info, lty=2, col='blue')
		plot(p)
		dev.off()		
	}
}










