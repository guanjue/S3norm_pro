### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
output_name = args[2]

sig = scan(signal_track_file)
thesh = 0

sig_non0 = sig[sig>thesh]
sig_mean = mean(sig_non0)
sig_var = var(sig_non0)
print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_var/sig_mean, digits=3)) ))

sig_prob = sig_mean / sig_var
if (sig_prob<0.1){
	print('use lower bound for nbp prob')
	sig_prob = 0.1
}

if (sig_prob>=0.9){
	print('use upper bound for nbp prob')
	sig_prob = 0.9
}

sig_size = sig_mean * sig_prob / (1-sig_prob)

nb_pval = apply(as.matrix(sig), MARGIN=1, function(x) pnbinom(x[1], sig_size, sig_prob, lower.tail=FALSE) )

### write output
write.table(nb_pval, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

