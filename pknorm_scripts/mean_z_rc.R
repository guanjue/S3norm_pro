library(metap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

cell_marker = args[1]
tail = args[2]
ctrl = args[3]
input_folder = args[4]
### extract filenames of the cell marker
file_list = list.files(input_folder, pattern=paste('^', cell_marker, '(.*)', tail, '$', sep='') )
print(file_list)
### read files of the cell marker
data_matrix = NULL
for (file in file_list){
	d = read.table(paste(input_folder, file, sep=''), header = F)
	data_matrix = cbind(data_matrix, d[,])
}

mean_sig = apply(data_matrix, 1, mean)

ctrl_sig = read.table(paste(input_folder, ctrl, sep=''), header = F)
mean_sig = (mean_sig+0.1) / (ctrl+0.1)

### write output
write.table(mean_sig, paste(cell_marker, '.mean_sig.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
mean_sig_log2 = log2(mean_sig)
z_sig = (mean_sig_log2-mean(mean_sig_log2))/sd(mean_sig_log2)
write.table(z_sig, paste(cell_marker, '.mean_sig_z.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### write FRiP & SNRs
pk_sig = mean_sig[p.adjust(pnorm(-z_sig),'fdr')<0.05]
bg_sig = mean_sig[p.adjust(pnorm(-z_sig),'fdr')>=0.05]
FRiP = sum(pk_sig) / sum(bg_sig)
SNR = mean(pk_sig) / mean(bg_sig)

write.table(c(FRiP, SNR), paste(cell_marker, '.mean_sig.frip_snr.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

