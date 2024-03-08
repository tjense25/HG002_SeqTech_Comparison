ml samtools bedtools

ICLR_BAM="./ICLR_HG002/grch38_vcf/Unmarked_Sample_611773165"
PB_BAM="./output/alignment/PacBio_Revio.GRCh38.sorted.bam"
ONT_BAM="./output/alignment/ONT_duplex.GRCh38.sorted.bam"
samtools depth -a -b ./truth_sets/GRCh38_mrg_full_exon.sorted.bed -Q 10 -J  $ICLR_BAM $PB_BAM $ONT_BAM > exon_depth.mapQ_10.txt 
#samtools depth -a -b ./truth_sets/GRCh38_mrg_full_gene.sorted.bed -Q 10 -J $ICLR_BAM $PB_BAM $ONT_BAM > gene_depth.mapQ_10.txt 

awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' exon_depth.mapQ_10.txt | bedtools intersect -loj -a - -b truth_sets/GRCh38_mrg_full_gene.sorted.bed | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' > ./output/exon_depth.mapQ_10.gene_annotated.txt 
#awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' gene_depth.mapQ_10.txt | bedtools intersect -loj -a - -b truth_sets/GRCh38_mrg_full_gene.sorted.bed | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' > ./output/gene_depth.mapQ_10.gene_annotated.txt 
