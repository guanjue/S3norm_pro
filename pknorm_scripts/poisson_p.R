### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
signal_folder = args[2]
input_track_file = args[3]
input_folder = args[4]
output_name = args[5]


get_poisson_pval = function(x){
	pval = ppois(x[1], lambda = x[2], lower.tail=FALSE)
	return(pval)
}

### read data
sig = read.table(paste(signal_folder, signal_track_file, sep=''), header = F)
input = read.table(paste(input_folder, input_track_file, sep=''), header = F)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
### get sig bg regions no bgs
thresh = 0

sig_od = sig[,1]
sig_mean = mean(sig_od)

### get local poisson mean
input_od = input[,1]
### use global mean if local poisson mean is smaller
input_od[input_od<sig_mean] = sig_mean

### get poisson p-value 1st round
sig_input = cbind(sig_od, input_od)
poisson_pval = apply(sig_input, MARGIN=1, function(x) get_poisson_pval(x) )

### get -log10(p-value)
poisson_pval[poisson_pval<=1e-324] = 1e-324
neglog10_poisson_pval = -log10(poisson_pval)

### write output
write.table(neglog10_poisson_pval, paste(output_name, '.poisson_pval.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

