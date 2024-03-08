library(tidyverse)
library(magrittr)


depths <- read_tsv("./output/exon_depth.mapQ_10.gene_annotated.txt", col_names=c("chrom", "pos", "ICLR", "PacBio", "ONT", "Gene")) 
percent_covered <- depths %>% pivot_longer(c(ICLR,ONT,PacBio), names_to="technology", values_to="depth") %>% group_by(Gene, technology) %>% summarize(percent_covered=sum(depth>10)/n()) %>% group_by(technology) %>% mutate(rank=rank(-percent_covered, ties.method = "first"))
percent_covered$technology %<>% factor(levels=c("ONT", "PacBio", "ICLR"))

summary(depths)
arrange(percent_covered, technology, rank) %>% 
 ggplot(aes(rank, percent_covered, color = technology)) + geom_line(size=2) +
    scale_color_manual(values=c("#0099CC", "#FF0066", "#FF9900")) + ylab("Percent of exonic bases covered > 10 depth") + xlab("Rank Order Gene") + theme_classic() + theme(panel.grid.major=element_line())
ggsave('./output/plots/CMRG_gene_percent_exon_coverage.depth_10.mapQ_10.all_genes.pdf')

not_fully_covered_rank <- percent_covered %>% filter(percent_covered < 1) %>% pull(rank) %>% min
arrange(percent_covered, technology, rank) %>% 
 ggplot(aes(rank, percent_covered, color  = technology)) + geom_line(size=2) + xlim(not_fully_covered_rank,NA) +
 scale_color_manual(values=c("#0099CC", "#FF0066", "#FF9900")) + ylab("Percent of exonic bases covered > 10 depth") + xlab("Rank Order Gene") + theme_classic() + theme(panel.grid.major=element_line())
ggsave('./output/plots/CMRG_gene_percent_exon_coverage.depth_10.mapQ_10.tail_genes.pdf')


