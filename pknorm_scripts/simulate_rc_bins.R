total_bin_num = 2e+5
pk_bins_num = 2e+4
dif_pk_bins_num = as.integer(pk_bins_num*0.2)
dif_fc = 4.0
random_distributed_rc = 0.95

all_bins = rep(0, total_bin_num)
### rep1 spe pk
all_bins[1:dif_pk_bins_num]=1
### rep2 spe pk
all_bins[(dif_pk_bins_num+1):(dif_pk_bins_num*2)]=2
### common pk
all_bins[(dif_pk_bins_num*2+1):pk_bins_num]=3

shape = m^2/v
rate = m/v














