args = commandArgs(trailingOnly=TRUE)

input_matrix = args[1]
output = args[2]
ylim_l = as.numeric(args[3])
ylim_h = as.numeric(args[4])
print(ylim_h)
data = read.table(input_matrix, header=TRUE, sep='\t')

if (is.na(ylim_l)){
	ylim_h = max(data[,-1])
	ylim_l = min(data[,-1])
}

pdf(output)
plot(data[,1], data[,4], pch=16, col='orchid1', ylim=c(ylim_l, ylim_h))
lines(data[,1], data[,4], lwd=1.5, col='orchid1')
points(data[,1], data[,2], pch=16, col='seagreen1')
lines(data[,1], data[,2], lwd=1.5, col='seagreen1')
points(data[,1], data[,3], pch=16, col='dodgerblue')
lines(data[,1], data[,3], lwd=1.5, col='dodgerblue')
points(data[,1], data[,5], pch=16, col='darkorange')
lines(data[,1], data[,5], lwd=1.5, col='darkorange')
points(data[,1], data[,6], pch=16, col='black')
lines(data[,1], data[,6], lwd=1.5, col='black')
dev.off()


