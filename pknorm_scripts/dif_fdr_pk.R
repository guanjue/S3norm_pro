args = commandArgs(trailingOnly=TRUE)

input_bed = args[1]
output_bed = args[2]

data = read.table(input_bed, header=FALSE, sep='\t')

data_p = data[,4]
data_p_fdr = p.adjust(data[,4], 'fdr')

thresh_list = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)

for (p_lim in thresh_list){
	data_pk = data[data_p_fdr<p_lim,]
	write.table(data_pk, paste(output_bed, '.', toString(p_lim), '.bed', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}

