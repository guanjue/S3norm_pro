cd /Users/universe/Documents/2018_BG/pknorm_analysis/h3k4me3/tss/
cd /storage/home/gzx103/scratch/vision/all_final_data/5end_reads_count/tss
cat /storage/home/gzx103/group/projects/vision/rna/rnaHtseqCountsall.txt | awk -F '\t' -v OFS='\t' '{print $1, $2+$3, $4+$5, $6+$7, $8+$9, $10+$11, $12+$13, $14+$15, $16+$17, $18+$18, $19+$20, $21+$22, $23+$24}' > rnaHtseqCountsall_replicate_merge.txt

cat /storage/home/gzx103/group/projects/vision/rna/gencode.vM4.annotation.bed | awk -F '\t' -v OFS='\t' '{if ($5=="protein_coding") print $1,$2,$3,$4,$5,$6}' > gencode.vM4.annotation.pc.bed


python ~/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/tss_sig/vlookup.py -t gencode.vM4.annotation.pc.bed -m 4 -s B_SPL.manorm.bed.tss.bed -n 4 -o gencode.vM4.annotation.pc.sorted.bed -k F

python ~/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/tss_sig/vlookup.py -t rnaHtseqCountsall_replicate_merge.txt -m 1 -s B_SPL.manorm.bed.tss.bed -n 4 -o rnaHtseqCountsall_replicate_merge.pcsorted.txt -k F


#gene	CFUE	CFUMk	CMP	ERY	GMP	iMK	LSK	MEP	MON	NEU	ER4	G1E

paste rnaHtseqCountsall_replicate_merge.pcsorted.txt gencode.vM4.annotation.pc.sorted.bed | awk -F '\t' -v OFS='\t' '{print $1, $2/($16-$15)*1000,  $4/($16-$15)*1000, $5/($16-$15)*1000, $6/($16-$15)*1000, $7/($16-$15)*1000, $8/($16-$15)*1000, $9/($16-$15)*1000, $10/($16-$15)*1000, $11/($16-$15)*1000, $12/($16-$15)*1000, $13/($16-$15)*1000}' > rna_rpk.pcsorted.txt

paste rnaHtseqCountsall_replicate_merge.pcsorted.txt gencode.vM4.annotation.pc.sorted.bed | awk -F '\t' -v OFS='\t' '{print $1, $2,  $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, ($16-$15)}' > rna_rc.pcsorted.txt

paste rnaHtseqCountsall_replicate_merge.pcsorted.txt gencode.vM4.annotation.pc.sorted.bed | awk -F '\t' -v OFS='\t' '{print $1, $2/($16-$15)*1000,  $3/($16-$15)*1000,  $4/($16-$15)*1000, $5/($16-$15)*1000, $6/($16-$15)*1000, $7/($16-$15)*1000, $8/($16-$15)*1000, $9/($16-$15)*1000, $10/($16-$15)*1000, $11/($16-$15)*1000, $12/($16-$15)*1000, $13/($16-$15)*1000}' > rna_rpk.pcsorted.all12.txt


paste CFU_E_ad.raw.bed.tss.txt CMP.raw.bed.tss.txt ERY_ad.raw.bed.tss.txt GMP.raw.bed.tss.txt MK_imm_ad.raw.bed.tss.txt LSK_BM.raw.bed.tss.txt MEP.raw.bed.tss.txt MONO_BM.raw.bed.tss.txt NEU.raw.bed.tss.txt ER4.raw.bed.tss.txt G1E.raw.bed.tss.txt > tss_h3k4me3.pcsorted.raw.txt
paste CFU_E_ad.raw.bed.1000kb.tss.txt CMP.raw.bed.1000kb.tss.txt ERY_ad.raw.bed.1000kb.tss.txt GMP.raw.bed.1000kb.tss.txt MK_imm_ad.raw.bed.1000kb.tss.txt LSK_BM.raw.bed.1000kb.tss.txt MEP.raw.bed.1000kb.tss.txt MONO_BM.raw.bed.1000kb.tss.txt NEU.raw.bed.1000kb.tss.txt ER4.raw.bed.1000kb.tss.txt G1E.raw.bed.1000kb.tss.txt > tss_h3k4me3.pcsorted.raw.1000kb.txt

paste CFU_E_ad.trnorm.bed.tss.txt CMP.trnorm.bed.tss.txt ERY_ad.trnorm.bed.tss.txt GMP.trnorm.bed.tss.txt MK_imm_ad.trnorm.bed.tss.txt LSK_BM.trnorm.bed.tss.txt MEP.trnorm.bed.tss.txt MONO_BM.trnorm.bed.tss.txt NEU.trnorm.bed.tss.txt ER4.trnorm.bed.tss.txt G1E.trnorm.bed.tss.txt > tss_h3k4me3.pcsorted.trnorm.txt
paste CFU_E_ad.trnorm.bed.1000kb.tss.txt CMP.trnorm.bed.1000kb.tss.txt ERY_ad.trnorm.bed.1000kb.tss.txt GMP.trnorm.bed.1000kb.tss.txt MK_imm_ad.trnorm.bed.1000kb.tss.txt LSK_BM.trnorm.bed.1000kb.tss.txt MEP.trnorm.bed.1000kb.tss.txt MONO_BM.trnorm.bed.1000kb.tss.txt NEU.trnorm.bed.1000kb.tss.txt ER4.trnorm.bed.1000kb.tss.txt G1E.trnorm.bed.1000kb.tss.txt > tss_h3k4me3.pcsorted.trnorm.1000kb.txt

paste CFU_E_ad.pknorm.bed.tss.txt CMP.pknorm.bed.tss.txt ERY_ad.pknorm.bed.tss.txt GMP.pknorm.bed.tss.txt MK_imm_ad.pknorm.bed.tss.txt LSK_BM.pknorm.bed.tss.txt MEP.pknorm.bed.tss.txt MONO_BM.pknorm.bed.tss.txt NEU.pknorm.bed.tss.txt ER4.pknorm.bed.tss.txt G1E.pknorm.bed.tss.txt > tss_h3k4me3.pcsorted.pknorm.txt
paste CFU_E_ad.pknorm.bed.1000kb.tss.txt CMP.pknorm.bed.1000kb.tss.txt ERY_ad.pknorm.bed.1000kb.tss.txt GMP.pknorm.bed.1000kb.tss.txt MK_imm_ad.pknorm.bed.1000kb.tss.txt LSK_BM.pknorm.bed.1000kb.tss.txt MEP.pknorm.bed.1000kb.tss.txt MONO_BM.pknorm.bed.1000kb.tss.txt NEU.pknorm.bed.1000kb.tss.txt ER4.pknorm.bed.1000kb.tss.txt G1E.pknorm.bed.1000kb.tss.txt > tss_h3k4me3.pcsorted.pknorm.1000kb.txt

paste CFU_E_ad.qtnorm.bed.tss.txt CMP.qtnorm.bed.tss.txt ERY_ad.qtnorm.bed.tss.txt GMP.qtnorm.bed.tss.txt MK_imm_ad.qtnorm.bed.tss.txt LSK_BM.qtnorm.bed.tss.txt MEP.qtnorm.bed.tss.txt MONO_BM.qtnorm.bed.tss.txt NEU.qtnorm.bed.tss.txt ER4.qtnorm.bed.tss.txt G1E.qtnorm.bed.tss.txt > tss_h3k4me3.pcsorted.qtnorm.txt
paste CFU_E_ad.qtnorm.bed.1000kb.tss.txt CMP.qtnorm.bed.1000kb.tss.txt ERY_ad.qtnorm.bed.1000kb.tss.txt GMP.qtnorm.bed.1000kb.tss.txt MK_imm_ad.qtnorm.bed.1000kb.tss.txt LSK_BM.qtnorm.bed.1000kb.tss.txt MEP.qtnorm.bed.1000kb.tss.txt MONO_BM.qtnorm.bed.1000kb.tss.txt NEU.qtnorm.bed.1000kb.tss.txt ER4.qtnorm.bed.1000kb.tss.txt G1E.qtnorm.bed.1000kb.tss.txt > tss_h3k4me3.pcsorted.qtnorm.1000kb.txt

paste CFU_E_ad.manorm.bed.tss.txt CMP.manorm.bed.tss.txt ERY_ad.manorm.bed.tss.txt GMP.manorm.bed.tss.txt MK_imm_ad.manorm.bed.tss.txt LSK_BM.manorm.bed.tss.txt MEP.manorm.bed.tss.txt MONO_BM.manorm.bed.tss.txt NEU.manorm.bed.tss.txt ER4.manorm.bed.tss.txt G1E.manorm.bed.tss.txt > tss_h3k4me3.pcsorted.manorm.txt
paste CFU_E_ad.manorm.bed.1000kb.tss.txt CMP.manorm.bed.1000kb.tss.txt ERY_ad.manorm.bed.1000kb.tss.txt GMP.manorm.bed.1000kb.tss.txt MK_imm_ad.manorm.bed.1000kb.tss.txt LSK_BM.manorm.bed.1000kb.tss.txt MEP.manorm.bed.1000kb.tss.txt MONO_BM.manorm.bed.1000kb.tss.txt NEU.manorm.bed.1000kb.tss.txt ER4.manorm.bed.1000kb.tss.txt G1E.manorm.bed.1000kb.tss.txt > tss_h3k4me3.pcsorted.manorm.1000kb.txt

