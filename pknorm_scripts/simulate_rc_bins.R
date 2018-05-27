total_bin_num = 200000
pk_bins_num = 20000
dif_pk_bins_num = 4000

all_bins = rep(0, total_bin_num)
### rep1 spe pk
all_bins[1:dif_pk_bins_num]=1
### rep2 spe pk
all_bins[(dif_pk_bins_num+1):(dif_pk_bins_num*2)]=2
### common pk
all_bins[(dif_pk_bins_num*2+1):pk_bins_num]=3


