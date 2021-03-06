ct=T_CD8_SPL
### get bedgraph file
tail -n+2 T_CD8_SPL.fisher_p.txt.all5info.txt.pkid.txt | awk -F '\t' -v OFS='\t' '{ if ($4==1) print $1,$2,$3, $10 }' | sort -k1,1 -k2,2n > ERY_fl.fisher_p.txt.all5info.txt.raw_pkid.bedgraph
cp ERY_fl.fisher_p.txt.all5info.txt.raw_pkid.bedgraph ERY_fl.fisher_p.txt.all5info.txt.qtnorm_pkid.bedgraph
cp ERY_fl.fisher_p.txt.all5info.txt.qtnorm_pkid.bedgraph ERY_fl.fisher_p.txt.all5info.txt.trnorm_pkid.bedgraph
cp ERY_fl.fisher_p.txt.all5info.txt.qtnorm_pkid.bedgraph ERY_fl.fisher_p.txt.all5info.txt.manorm_pkid.bedgraph
cp ERY_fl.fisher_p.txt.all5info.txt.qtnorm_pkid.bedgraph ERY_fl.fisher_p.txt.all5info.txt.pknorm_pkid.bedgraph


for ct in $(cat ct_list.txt)
do

ct=T_CD8_SPL

tail -n+2 $ct'.fisher_p.txt.all5info.txt.pkid.txt' | awk -F '\t' -v OFS='\t' '{ if ($5==1) print $1,$2,$3, $11 }' | sort -k1,1 -k2,2n > $ct'.fisher_p.txt.all5info.txt.raw_pkid.bedgraph'
tail -n+2 $ct'.fisher_p.txt.all5info.txt.pkid.txt' | awk -F '\t' -v OFS='\t' '{ if ($6==1) print $1,$2,$3, $12 }' | sort -k1,1 -k2,2n > $ct'.fisher_p.txt.all5info.txt.trnorm_pkid.bedgraph'
tail -n+2 $ct'.fisher_p.txt.all5info.txt.pkid.txt' | awk -F '\t' -v OFS='\t' '{ if ($7==1) print $1,$2,$3, $13 }' | sort -k1,1 -k2,2n > $ct'.fisher_p.txt.all5info.txt.manorm_pkid.bedgraph'
tail -n+2 $ct'.fisher_p.txt.all5info.txt.pkid.txt' | awk -F '\t' -v OFS='\t' '{ if ($8==1) print $1,$2,$3, $14 }' | sort -k1,1 -k2,2n > $ct'.fisher_p.txt.all5info.txt.qtnorm_pkid.bedgraph'
tail -n+2 $ct'.fisher_p.txt.all5info.txt.pkid.txt' | awk -F '\t' -v OFS='\t' '{ if ($9==1) print $1,$2,$3, $15 }' | sort -k1,1 -k2,2n > $ct'.fisher_p.txt.all5info.txt.pknorm_pkid.bedgraph'
### get different pks
bedtools intersect -a $ct'.fisher_p.txt.all5info.txt.qtnorm_pkid.bedgraph' -b $ct'.fisher_p.txt.all5info.txt.pknorm_pkid.bedgraph' -v > $ct'.qt_nopk.bed'
bedtools intersect -a $ct'.fisher_p.txt.all5info.txt.trnorm_pkid.bedgraph' -b $ct'.fisher_p.txt.all5info.txt.pknorm_pkid.bedgraph' -v > $ct'.tr_nopk.bed'
bedtools intersect -a $ct'.fisher_p.txt.all5info.txt.raw_pkid.bedgraph' -b $ct'.fisher_p.txt.all5info.txt.pknorm_pkid.bedgraph' -v > $ct'.raw_nopk.bed'
bedtools intersect -a $ct'.fisher_p.txt.all5info.txt.pknorm_pkid.bedgraph' -b $ct'.fisher_p.txt.all5info.txt.manorm_pkid.bedgraph' -v > $ct'.pk_noma.bed'

tp_list=(qt_nopk tr_nopk raw_nopk pk_noma)
for i in {0..3}
do
	echo ${tp_list[i]}
	### getfasta
	bedtools getfasta -fi ~/group/genome/mm10/mm10_no_alt_analysis_set_ENCODE.fasta -bed $ct'.'${tp_list[i]}'.bed' > $ct'.'${tp_list[i]}'.fa'
	### get fimo result
	time ~/group/software/meme/bin/fimo --text ctcf_motif.meme $ct'.'${tp_list[i]}'.fa' > $ct'.'${tp_list[i]}'.fimo.txt'
	### get fimo bed
	tail -n+2 $ct'.'${tp_list[i]}'.fimo.txt' | cut -f3 | awk -F ':' -v OFS='\t' '{print $1,$2}' | awk -F '-' -v OFS='\t' '{print $1,$2}' > $ct'.'${tp_list[i]}'.fimo.bed'
	### get intersect fimo motif bed
	bedtools intersect -a $ct'.'${tp_list[i]}'.bed' -b $ct'.'${tp_list[i]}'.fimo.bed' -wa > $ct'.'${tp_list[i]}'.fimomotif.bed'
	bedtools merge -i $ct'.'${tp_list[i]}'.fimomotif.bed' > $ct'.'${tp_list[i]}'.merge.fimomotif.bed'
	wc -l $ct'.'${tp_list[i]}'.merge.fimomotif.bed'
	bedtools merge -i $ct'.'${tp_list[i]}'.bed' > $ct'.'${tp_list[i]}'.merge.bed'
	wc -l $ct'.'${tp_list[i]}'.merge.bed'
	### get NOT intersect fimo motif bed
	bedtools intersect -a $ct'.'${tp_list[i]}'.bed' -b $ct'.'${tp_list[i]}'.fimo.bed' -v > $ct'.'${tp_list[i]}'.NOfimomotif.bed'
	bedtools merge -i $ct'.'${tp_list[i]}'.NOfimomotif.bed' > $ct'.'${tp_list[i]}'.merge.NOfimomotif.bed'
done

done






mt_list=(raw_pkid trnorm_pkid manorm_pkid qtnorm_pkid pknorm_pkid)
for i in {0..4}
do
	echo ${mt_list[i]}
	### getfasta
	bedtools getfasta -fi ~/group/genome/mm10/mm10_no_alt_analysis_set_ENCODE.fasta -bed $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.bedgraph' > $ct'.'${mt_list[i]}'.fa'
	### get fimo result
	time ~/group/software/meme/bin/fimo --text ctcf_motif.meme $ct'.'${mt_list[i]}'.fa' > $ct'.'${mt_list[i]}'.fimo.txt'
	### get fimo bed
	tail -n+2 $ct'.'${mt_list[i]}'.fimo.txt' | cut -f3 | awk -F ':' -v OFS='\t' '{print $1,$2}' | awk -F '-' -v OFS='\t' '{print $1,$2}' > $ct'.'${mt_list[i]}'.fimo.bed'
	### get intersect fimo motif bed
	bedtools intersect -a $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.bedgraph' -b $ct'.'${mt_list[i]}'.fimo.bed' -wa > $ct'.'${mt_list[i]}'.fimomotif.bed'
	bedtools merge -i $ct'.'${mt_list[i]}'.fimomotif.bed' > $ct'.'${mt_list[i]}'.merge.fimomotif.bed'
	wc -l $ct'.'${mt_list[i]}'.merge.fimomotif.bed'
	bedtools merge -i $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.bedgraph' > $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.merge.bedgraph'
	wc -l $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.merge.bedgraph'
	### get peak bigwig
	#~/group/software/ucsc/bedGraphToBigWig $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.merge.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.fisher_p.txt.all5info.txt.'${mt_list[i]}'.bw'
done

time bedtools map -a T_CD8_SPL.qt_nopk.merge.NOfimomotif.bed -b T_CD8_SPL.fisher_p.txt.all5info.txt.qtnorm_pkid.all.bedgraph -c 4 -o mean > T_CD8_SPL.qt_nopk.merge.NOfimomotif.qtmeansig.bed
time bedtools map -a T_CD8_SPL.qt_nopk.merge.NOfimomotif.bed -b T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_pkid.all.bedgraph -c 4 -o mean > T_CD8_SPL.qt_nopk.merge.NOfimomotif.pkmeansig.bed

paste T_CD8_SPL.qt_nopk.merge.NOfimomotif.qtmeansig.bed T_CD8_SPL.qt_nopk.merge.NOfimomotif.pkmeansig.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,($4+10.0)/($8+10.0)}' | sort -k4,4nr | awk -F '\t' -v OFS='\t' '{print $1":"$2"-"$3,$4}' > T_CD8_SPL.qt_nopk.merge.NOfimomotif.qtpkmeansig.bed

time bedtools map -a T_CD8_SPL.pk_noma.merge.fimomotif.bed -b T_CD8_SPL.fisher_p.txt.all5info.txt.manorm_pkid.all.bedgraph -c 4 -o mean > T_CD8_SPL.pk_noma.merge.fimomotif.mameansig.bed
time bedtools map -a T_CD8_SPL.pk_noma.merge.fimomotif.bed -b T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_pkid.all.bedgraph -c 4 -o mean > T_CD8_SPL.pk_noma.merge.fimomotif.pkmeansig.bed

paste T_CD8_SPL.pk_noma.merge.fimomotif.pkmeansig.bed T_CD8_SPL.pk_noma.merge.fimomotif.mameansig.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,($4+10.0)/($8+10.0)}' | sort -k4,4nr | awk -F '\t' -v OFS='\t' '{print $1":"$2"-"$3,$4}' > T_CD8_SPL.pk_noma.merge.fimomotif.pkmameansig.bed





time computeMatrix scale-regions -S T_CD8_SPL.manorm.sort.bw T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_all.sig.bw -R T_CD8_SPL.pk_noma.merge.fimomotif.bed -a 0 -b 0 --binSize 200 -out T_CD8_SPL.pk_noma.merge.fimomotif.pkmameansig.mat.gz


cat T_CD8_SPL.fisher_p.txt.all5info.txt.raw_pkid.merge.bedgraph | awk '{print $1,$2,$3,1}' > T_CD8_SPL.fisher_p.txt.all5info.txt.raw_pkid.merge.1.bedgraph

~/group/software/ucsc/bedGraphToBigWig T_CD8_SPL.fisher_p.txt.all5info.txt.raw_pkid.merge.1.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes T_CD8_SPL.raw_pkid.sort.bw

~/group/software/ucsc/bedGraphToBigWig T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_pkid.merge.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_pkid.merge.bw
~/group/software/ucsc/bedGraphToBigWig T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_pkid.merge.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes T_CD8_SPL.fisher_p.txt.all5info.txt.pknorm_pkid.merge.bw








