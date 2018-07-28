args = commandArgs(trailingOnly=TRUE)
input_sig1 = args[1]
output = args[2]

small_num = 1
random_sample_num = 1000000
upperlim = 1000
lowerlim = 0
ct = unlist(strsplit(input_sig1, "[.]"))[1]

sig1 = scan(input_sig1)
sig1_log2 = log2(sig1+small_num)
sig1_z = (sig1_log2 - mean(sig1_log2))/sd(sig1_log2)
#sig1_z = (sig1 - mean(sig1))/sd(sig1)

write.table(sig1_z, paste(output,".znorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)





