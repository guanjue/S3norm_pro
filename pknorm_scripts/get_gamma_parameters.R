#get_gamma_parameters.R

file_list = list.files(pattern = "\\.bamtobed5endintersect.signal$")

d_all = c()
for (filename in file_list){
	print(filename)
	d_tmp = scan(filename)
	d_all = cbind(d_all, d_tmp)
}

ctcf_col_means = colMeans(d_all)
ctcf_col_means_mean = mean(ctcf_col_means)

d_all_norm = t(apply(d_all, 1, function(x) x/ctcf_col_means*ctcf_col_means_mean))
ctcf_col_vars = apply(d_all_norm, 1, var)

d_all_norm_fc = c()
for (i in c(1:dim(d_all_norm)[2])){
	for (j in c(1:dim(d_all_norm)[2])){
		print(paste(toString(i), '_', toString(j), sep=''))
		if (i < j){
			d_all_norm_fc_tmp = (d_all_norm[,i]+0.01) / (d_all_norm[,j]+0.01)
			d_all_norm_fc = cbind(d_all_norm_fc, d_all_norm_fc_tmp)
		}
	}
}

###### fold change mean & variance
d_all_norm_fc_log2 = sort(as.vector(log2(d_all_norm_fc)))
top_rank = length(d_all_norm_fc_log2) * 0.025

ctcf_fc_log2_up_means = mean(head(d_all_norm_fc_log2, top_rank))
ctcf_fc_log2_up_vars = var(head(d_all_norm_fc_log2, top_rank))
ctcf_fc_log2_down_means = mean(tail(d_all_norm_fc_log2, top_rank))
ctcf_fc_log2_down_vars = var(tail(d_all_norm_fc_log2, top_rank))
ctcf_fc_log2_nochange_means = mean(d_all_norm_fc_log2[(top_rank+1):(length(d_all_norm_fc_log2)-top_rank)])
ctcf_fc_log2_nochange_vars = var(d_all_norm_fc_log2[(top_rank+1):(length(d_all_norm_fc_log2)-top_rank)])

fc_list = cbind(c(ctcf_fc_log2_nochange_means,ctcf_fc_log2_nochange_vars), c(ctcf_fc_log2_up_means,ctcf_fc_log2_up_vars), c(ctcf_fc_log2_down_means,ctcf_fc_log2_down_vars))

fc_list = cbind(c(0,ctcf_fc_log2_nochange_vars), c(2.0,ctcf_fc_log2_up_vars), c(-2.0,ctcf_fc_log2_down_vars))

fc_list[2,]=sqrt(fc_list[2,])

###### mean mean & variance
ctcf_all_means = mean(d_all_norm)
ctcf_all_vars = var(as.vector(d_all_norm))

###### gamma mean shape & rate
ctcf_shape = ctcf_all_means^2 / ctcf_all_vars
ctcf_rate = ctcf_all_means / ctcf_all_vars

x = seq(0, 100, length=1000)
hx = dgamma(x, shape = ctcf_shape, rate = ctcf_rate)

png('test_gamma_dist.png')
plot(x, hx, type="l", lty=2, xlab="x value",
  ylab="Density", main="Comparison of t Distributions")
dev.off()


samplesize = 1e+5
set.seed(2018)

porp <- c(0.04,0.95,0.01);  # nonD, upperD, lowerD

dsimu <- matrix(0,nrow=samplesize,ncol=3);
temp1 <- runif(samplesize,0,1);        # determine DEGs
dsimu[temp1<porp[1],1] <- 1;           # 1 for upper
dsimu[temp1>(porp[1]+porp[2]),1] <- 2; # 2 for lower

dispersion = 0.1
fc_used=c()
mu_used=c()

for (simui in 1:samplesize){
	if (simui%%(100000)==0){
		print(simui)
	}
	rdq_mu  <- runif(1,0,1);      # random quantail for overall mean
	rdq_lfc <- runif(1,0,1);      # random quantail for log fold change
	# for seq
	sp.s.NBmu <- qgamma(rdq_mu,shape=ctcf_shape,rate=ctcf_rate);
	mu_used[simui] = sp.s.NBmu
	fc_used[simui] = 2^(qnorm(rdq_lfc,mean=fc_list[1,dsimu[simui,1]+1],sd=fc_list[2,dsimu[simui,1]+1])/2)
	sp.s.NBmu_up <- sp.s.NBmu*2^(qnorm(rdq_lfc,mean=fc_list[1,dsimu[simui,1]+1],sd=fc_list[2,dsimu[simui,1]+1])/2)
	sp.s.NBmu_down <- sp.s.NBmu*2^(-qnorm(rdq_lfc,mean=fc_list[1,dsimu[simui,1]+1],sd=fc_list[2,dsimu[simui,1]+1])/2);
	sp.s.NBsize_up <- sp.s.NBmu_up/9#1/rgamma(1,shape=p.seq$dgam[1],rate=p.seq$dgam[2]);   # size for gamma mean dispersion parameter, little different from the one in paper (1/r)
	sp.s.NBsize_down <- sp.s.NBmu_down/9
	dsimu[simui,2] <- rnbinom(1,mu=sp.s.NBmu_up,size=0.1);
	dsimu[simui,3] <- rnbinom(1,mu=sp.s.NBmu_down,size=0.1);
}

sample_id = sample(dim(dsimu)[1], 100000)

#png('test_sim.png')
plot(dsimu[sample_id,2]+0.01, dsimu[sample_id,3]+0.01, log='xy', xlim=c(0.01, 1e+5), ylim=c(0.01, 1e+5))
abline(0,1,col='red')
#dev.off()


library(LSD)
png('test_sim_heat.png')
heatscatter(dsimu[,2], dsimu[,3], pch=16, log='xy', xlim=c(0.01, 1e+5), ylim=c(0.01, 1e+5), cex=0.6)
abline(0,1,col='red')
dev.off()







