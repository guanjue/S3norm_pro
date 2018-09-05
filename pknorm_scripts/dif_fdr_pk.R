args = commandArgs(trailingOnly=TRUE)

input_bed = args[1]
output_bed = args[2]

data = read.table(input_bed, header=FALSE, sep='\t')

data_p = 10^(-data[,4])
data_p_fdr = data_p

#data_p_fdr = p.adjust(data_p, 'fdr')

thresh_list = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 5e-10, 1e-10)

for (p_lim in thresh_list){
	print(p_lim)
	data_pk = data[data_p_fdr<p_lim,]
	write.table(data_pk, paste(output_bed, '.', toString(p_lim), '.bed', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}

