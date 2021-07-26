chrom3d-fullgenome
#1

mkdir -p lpshic/bedpe/intra
mkdir -p lpshic/bedpe/inter
mkdir -p lpshic/matrix
mkdir  -p lpshic/tads

#saline same

#2
scp npc_1038000000_merge_ko* yuezhu@192.168.79.84:/home/yuezhu/chrom3d/npc_hic_merge/ko
#saline same

#3  Convert intrachromosomal Hi-C to BEDPE an matrix format
# BEDPE:
awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }' /home/yuezhu/chrom3d/npc_hic_merge/wt/npc_1038000000_merge_wt_50000_abs.bed /home/yuezhu/chrom3d/npc_hic_merge/wt/npc_1038000000_merge_wt_50000.matrix | awk '$4==$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3>"salinehic/bedpe/intra/"$4}'

awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }' /home/yuezhu/chrom3d/npc_hic_merge/ko/npc_1038000000_merge_ko_50000_abs.bed /home/yuezhu/chrom3d/npc_hic_merge/ko/npc_1038000000_merge_ko_50000.matrix | awk '$4==$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3>"lpshic/bedpe/intra/"$4}'



# Matrix format:
for chr in salinehic/bedpe/intra/*
do
chrname=$(basename $chr)
cut -f 2,5,7 $chr > salinehic/matrix/$chrname
done

for chr in lpshic/bedpe/intra/*
do
chrname=$(basename $chr)
cut -f 2,5,7 $chr > lpshic/matrix/$chrname
done

#4 Convert interchromosomal Hi-C to BEDPE
awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }'  /home/yuezhu/chrom3d/npc_hic_merge/wt/npc_1038000000_merge_wt_1000000_abs.bed /home/yuezhu/chrom3d/npc_hic_merge/wt/npc_1038000000_merge_wt_1000000.matrix  | awk '$4<$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3>"salinehic/bedpe/inter/"$4"_"$7}'


awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }'  /home/yuezhu/chrom3d/npc_hic_merge/ko/npc_1038000000_merge_ko_1000000_abs.bed /home/yuezhu/chrom3d/npc_hic_merge/ko/npc_1038000000_merge_ko_1000000.matrix  | awk '$4<$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3>"lpshic/bedpe/inter/"$4"_"$7}'


#5.Convert called TADs into a segmented genome to define Chrom3D beads

bedtools complement -i lpshic/tads/npc_50000_allko.h5.bed -g /home/yuezhu/yard/HiC-Pro-3.0.0/annotation/chrom_mm10.sizes | cat - lpshic/tads/npc_50000_allko.h5.bed | bedtools sort -g /home/yuezhu/yard/HiC-Pro-3.0.0/annotation/chrom_mm10.sizes | awk '$1!="chrY"' | awk '$1!="chrM"' > npc_50000_allko.h5_beads.bed


6.Map intra-chromosomal interactions from Hi-C to the beads defined in the previous step and aggregate the contacts between these beads

cat salinehic/bedpe/intra/chr* | awk '{printf("%s\t%s\t%s\n",$1,$2,$2+1)}' | bedtools intersect -wao -a stdin -b npc_50000_allwt.h5_beads.bed | cut -f 4,5,6 > left.tmp
cat salinehic/bedpe/intra/chr* | awk '{printf("%s\t%s\t%s\t%s\n",$4,$5,$5+1,$7)}' | bedtools intersect -wao -a stdin -b npc_50000_allwt.h5_beads.bed | awk '{printf("%s\t%s\t%s\t%s\n",$5,$6,$7,$4)}' > right.tmp

paste left.tmp right.tmp | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6] += $7} END{for (i in a) print i"\t"a[i]}' |  awk '$2!=$5' | sort -k 2n,2n > npc_50000_allwt.h5_bead_interactions.intra.bedpe
rm left.tmp right.tmp


7 Remove interactions between beads overlapping centromeres

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | bedtools pairtobed -a D0_bead_interactions.intra.bedpe -b stdin -type neither > D0_bead_interactions.intra.nocen.bedpe

找不到mm10着丝粒bedfile,明天再找找，0712做到这儿了
0712：找不到，第七步不做了




8. Identifying statistically significant inter-bead interactions within chromosomes, using the Non-central Hypergeometric distribution (NCHG)

processing_scripts/NCHG_hic/NCHG -m 50000 -p npc_50000_allko.h5_bead_interactions.intra.bedpe > npc_50000_allko.h5_bead_interactions.intra.NCHG.out

# Correcting for multiple-testing using FDR:
conda activate py3
python processing_scripts/NCHG_fdr_oddratio_calc.py npc_50000_allko.h5_bead_interactions.intra.NCHG.out fdr_bh 2 0.01 > npc_50000_allko.h5_bead_interactions.intra.NCHG.sig

processing_scripts/NCHG_hic/NCHG -m 50000 -p npc_50000_allwt.h5_bead_interactions.intra.bedpe > npc_50000_allwt.h5_bead_interactions.intra.NCHG.out



9. Identifying statistically significant inter-bead interactions between chromosomes

# Create a BEDPE file containing interchromosomal interactions where 'blacklisted' regions are removed:

cat salinehic/bedpe/inter/chr* | bedtools pairtobed -type neither -a stdin -b mm10.blacklist.bed | python processing_scripts/cap_chr_end.py chrom_mm10.sizes > npc_50000_allwt.h5_bead_interactions.inter.noblist.bedpe


# Run NCHG, like above, but on the interchromosomal interactions: 
processing_scripts/NCHG_hic/NCHG -i -p npc_50000_allwt.h5_bead_interactions.inter.noblist.bedpe > npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.out


python processing_scripts/NCHG_fdr_oddratio_calc.py npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.out fdr_bh 2 0.01 > npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.sig







10 Mapping all (intra- and interchromosomal interactions) to the TADs/beads


awk '{printf("%s\t%s\t%s\n",$1,($2+$3)/2,1+($2+$3)/2)}' npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.sig >wtleftstdin

stdin remove e
filechange to STDIN.BED 

bedtools intersect -wao -a wtleftstdin.bed -b npc_50000_allwt.h5_beads.bed | cut -f 4,5,6 > left.tmp



awk '{printf("%s\t%s\t%s\t%s\n",$4,($5+$6)/2,1+($5+$6)/2,$7)}' npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.sig >wtrightstdin
| bedtools intersect -wao -a wtrightstdin.bed -b npc_50000_allwt.h5_beads.bed | awk '{printf("%s\t%s\t%s\t%s\n",$5,$6,$7,$4)}' > right.tmp

paste left.tmp right.tmp | sort -u -k 2n,2n > npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.tadwise.sig
rm left.tmp right.tmp


cat npc_50000_allwt.h5_bead_interactions.intra.NCHG.sig npc_50000_allwt.h5_bead_interactions.inter.noblist.NCHG.tadwise.sig | awk '$2!=-1 && $5!=-1' > npc_50000_allwt.h5_interactions.all.sig


wt&ko done!


11 Generate the Chrom3D input file in GTrack format, specifying the 3D model setup

python processing_scripts/makeGtrack.py npc_50000_allko.h5_interactions.all.sig npc_50000_allko.h5_beads.bed > npc_50000_allko.h5_interactions.gtrack



python processing_scripts/makeGtrack.py npc_50000_allwt.h5_interactions.all.sig npc_50000_allwt.h5_beads.bed > npc_50000_allwt.h5_interactions.gtrack


echo -e "##gtrack version: 1.0\n##track type: linked segments\n###seqid\tstart\tend\tid\tradius\tedges\n/" > head.gtrack

cat head.gtrack  npc_50000_allwt.h5_interactions.gtrack >npc_50000_allwt.h5_interactions.final.gtrack



13 Make the GTrack file define a diploid genome structure (optional)
python processing_scripts/make_diploid_gtrack.py npc_50000_allwt.h5_interactions.final3.gtrack > npc_50000_allwt.h5_interactions.final3.diploid.gtrack


sample:
##gtrack version: 1.0
##track type: linked segments
###seqid	start	end	id	radius	edges
chr1_A	0	3700000	chr1_A:0-3700000	1	chr1_A:87200000-88300000;chr10_A:5850000-6800000
chr1_B	0	3700000	chr1_B:0-3700000	1	chr1_B:87200000-88300000;chr10_B:5850000-6800000
chr1_A	3700000	4350000	chr1_A:3700000-4350000	1	chr1_A:87200000-88300000
chr1_B	3700000	4350000	chr1_B:3700000-4350000	1	chr1_B:87200000-88300000
chr1_A	4350000	4800000	chr1_A:4350000-4800000	1	chr1_A:87200000-88300000;chr12_A:17850000-20150000
chr1_B	4350000	4800000	chr1_B:4350000-4800000	1	chr1_B:87200000-88300000;chr12_B:17850000-20150000
chr1_A	4800000	5150000	chr1_A:4800000-5150000	1	chr1_A:87200000-88300000




14 Run Chrom3D based on the GTrack file (takes up to 20 hrs)
#! /bin/bash
#PBS -N chrome3dko
#PBS -l nodes=1
#PBS -walltime=1-3:00:00
#PBS -o chrom3dko.out
#PBS -e chrom3dko.err
cd /home/yuezhu/chrom3d/npc_hic_merge
Chrom3D -l 10000 -y 0.15 -r 5.0 -n 3000000 npc_50000_allko.h5_interactions.final3.diploid.gtrack > npc_50000_allko.h5_interactions.final3.diploid_model_full.cmm 2> npc_50000_allko.h5_interactions.final3.diploid.gtrack_model_full.err 


0716

color beads

unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/color_beads.py

ids sample:
awk '$6==1' D0_bead_interactions.lads.gtrack | cut -f 4 > lads.ids
head lads.ids 
chr18:700000-2650000
chr18:4000000-5050000
chr18:5050000-7200000
chr18:23700000-24900000
chr18:24900000-28650000
chr18:28650000-28750000
chr18:28750000-28850000
chr18:28850000-29050000
chr18:29050000-29500000
chr18:29500000-29700000


python3 processing_scripts/color_beads.py npc_50000_allko.h5_interactions.final3.diploid_model_full.cmm npc_ko_chr1.ids 255,0,0 OVERRIDE > npc_50000_allko.h5_interactions.final3.diploid_model_full_redchr1.cmm

python3 processing_scripts/color_beads.py npc_50000_allwt.h5_interactions.final3.diploid_model_full.cmm npc_wt_chr4.ids 255,0,0 OVERRIDE > npc_50000_allwt.h5_interactions.final3.diploid_model_full_redchr4.cmm

python3 processing_scripts/color_beads.py npc_50000_allko.h5_interactions.final3.diploid_model_full.cmm npc_ko_chr4.ids 255,0,0 OVERRIDE > npc_50000_allko.h5_interactions.final3.diploid_model_full_redchr4.cmm











