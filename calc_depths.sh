ml samtools bedtools

#samtools depth -a -b ./truth_sets/GRCh38_mrg_full_exon.sorted.bed -Q 5 ./output/alignment/ICLR.GRCh38.sorted.bam ./output/alignment/PacBio_Revio.GRCh38.sorted.bam ./output/alignment/ONT_duplex.GRCh38.sorted.bam > exon_depth.mapQ_5.txt &
#samtools depth -a -b ./truth_sets/GRCh38_mrg_full_gene.sorted.bed -Q 5 ./output/alignment/ICLR.GRCh38.sorted.bam ./output/alignment/PacBio_Revio.GRCh38.sorted.bam ./output/alignment/ONT_duplex.GRCh38.sorted.bam > gene_depth.mapQ_5.txt &
#wait

awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' exon_depth.mapQ_5.txt | bedtools intersect -loj -a - -b truth_sets/GRCh38_mrg_full_gene.sorted.bed | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' > ./output/exon_depth.mapQ_5.gene_annotated.txt &
awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' gene_depth.mapQ_5.txt | bedtools intersect -loj -a - -b truth_sets/GRCh38_mrg_full_gene.sorted.bed | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' > ./output/gene_depth.mapQ_5.gene_annotated.txt &
wait
