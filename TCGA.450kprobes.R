library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(tidyverse)

data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


methylation.results <- read.csv("../../Tables/Results/Methylation.Final.merged.table.wo.vus.with.mut.anot.csv")


str(annotation.table)

rownames(annotation.table)


annotation.df <- data.frame(rownames(annotation.table), annotation.table$chr, annotation.table$pos, 
                            annotation.table$UCSC_RefGene_Name, annotation.table$HMM_Island, annotation.table$Enhancer, 
                            annotation.table$Regulatory_Feature_Group)

merge.table <- inner_join(methylation.results, annotation.df, by = c("Gene.name"="rownames.annotation.table.") )

merge.table$START <- merge.table$annotation.table.pos
merge.table$END  <- merge.table$annotation.table.pos + 1

#bed file

bed.table <- merge.table %>% select(annotation.table.chr, START, END, Gene.name)
write.table (bed.table, "TCGA.450kprobes.bed", sep = "\t" , quote = F, row.names = F, col.names = F)


write.csv(merge.table, "Final.merged.table.methylation.annotated.with.mutation.annot.csv")


methylation.with.counts <- read.table("TCGA.450kprobes.ext.500.w.counts.MITF.REST.ARID2.FOSL2.bed", sep = "\t")

methylation.final <- inner_join(methylation.with.counts,merge.table, by = c("V4" ="Gene.name") )

methylation.final <- methylation.final %>%  mutate(methylation.status = ifelse(p.value_bh < 0.05 & Estimate_phenotype > 0.5, 'gain', "constant"))


ggplot(methylation.final, aes( x = methylation.status, y = log2(V7 + 1), fill = methylation.status )) +
  geom_boxplot()
