library(mixtools)
library(mclust)

set.seed(2018)
### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
output_name = args[2]

sig = scan(signal_track_file)
thesh = -1

used_id = sig>=thesh

sig_noLow_raw = sig[used_id]
sig_noLow_raw_addnorm = sig_noLow_raw+rnorm(length(sig_noLow_raw),1,0.1)
sig_noLow_raw_addnorm[sig_noLow_raw_addnorm<1]=1
sig_noLow = log(sig_noLow_raw_addnorm)

cluster_num = 3

if (length(sig_noLow)>100000){
	sample_num = 100000
	sample_id = sample(length(sig_noLow), sample_num)
	sig_noLow_s = sig_noLow[sample_id]
	sig_used = sig_noLow_s
}

sig_used = sig_noLow
mfit=Mclust(sig_used,G=cluster_num)


#mixmdl2 = normalmixEM(sig_used, k=cluster_num)

c1 = mfit$classification == 1
summary(sig_used[c1])
c2 = mfit$classification == 2
summary(sig_used[c2])
c3 = mfit$classification == 3
summary(sig_used[c3])

c_order = order(mfit$parameters$mean)
c_used = c(1:cluster_num)[c_order]
colors = c('red', 'green', 'blue')[c_order]

x=seq(0,max(sig_noLow),0.02)

png(paste(output_name, '.png', sep=''))

hist(sig_used, breaks=50, freq=FALSE)
s1 = (mfit$parameters$variance$sigmasq[1])^0.5
m1 = mfit$parameters$mean[1]
l1 = mfit$parameters$pro[1]
hx1 <- dnorm(x,m1,s1)
lines(x, hx1*l1, type='l', lty=1, col=colors[1])
s2 = (mfit$parameters$variance$sigmasq[2])^0.5
m2 = mfit$parameters$mean[2]
l2 = mfit$parameters$pro[2]
hx2 <- dnorm(x,m2,s2)
lines(x, hx2*l2, type='l', lty=1, col=colors[2])
s3 = (mfit$parameters$variance$sigmasq[3])^0.5
m3 = mfit$parameters$mean[3]
l3 = mfit$parameters$pro[3]
hx3 <- dnorm(x,m3,s3)
lines(x, hx3*l3, type='l', lty=1, col=colors[3])

dev.off()

c_order = order(mixmdl2$mu)
c_used = c(1:cluster_num)[order(mixmdl2$mu)]
#print(mixmdl2$mu)
#print(c_used)
pk2 = sig_used[apply(mixmdl2$posterior, 1, function(x) max(x)==x[c_used[cluster_num]])]
pk1 = sig_used[apply(mixmdl2$posterior, 1, function(x) max(x)==x[c_used[cluster_num-1]])]
pk0 = sig_used[apply(mixmdl2$posterior, 1, function(x) max(x)==x[c_used[cluster_num-2]])]

pk1_lim = min(exp(pk1))
pk2_lim = min(exp(pk2))
#print(pk1_lim)
#print(pk2_lim)

#print(sum(as.logical((sig<pk2_lim)*(sig>=pk1_lim))))

gmm_trinary = rep(3,length(sig))
gmm_trinary[sig>=pk2_lim] = 1
gmm_trinary[as.logical((sig<pk2_lim)*(sig>=pk1_lim))] = 0

### write output
write.table(gmm_trinary, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

