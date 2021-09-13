##
## File: plot-gnn.R
## Created: Fri Mar 12 20:11:07 2021
## $Id: $
#
## Copyright (C) 2021 by Per Unneberg
##
## Author: Per Unneberg
##
library(tidyverse)
library(ggthemes)

df.gnn <- read.csv(snakemake@input[["gnn"]])
samples <- snakemake@params[["samples"]]
df.samples <- read_tsv(snakemake@input[["samples"]])
df.pops <- read_tsv(snakemake@input[["populations"]])
df.meta.all <- merge(df.samples, df.pops)
df.meta <- df.meta.all[, c("SM", "population")]
colnames(df.meta) <- c("vcf.id", "taxa")
df.meta$taxa <- gsub("_[0-9]+$", "", df.meta$taxa)
title <- paste0("chromosome ", snakemake@wildcards[["chrom"]])

df.gnn.long <- df.gnn %>%
  mutate(haplo = as.factor(paste(Sample.node, Individual, sep = "_"))) %>%
  select(-Sample.node) %>%
  pivot_longer(names_to = "popGroup", values_to = "prob", cols = -c(haplo, Species, Individual)) %>%
  dplyr::rename(test.species = Species)

df.gnn.long <- df.gnn.long %>% left_join(df.meta, by = c("Individual" = "vcf.id"))

#create sorting
haplo.order <- df.gnn.long %>%
  filter(taxa == popGroup) %>%
  group_by(taxa) %>%
  arrange(-prob) %>% ungroup() %>% select(haplo) %>%
  mutate(order = 1:n())

df.gnn.long <- df.gnn.long %>% left_join(haplo.order, by = "haplo")

library(RColorBrewer)
nb.cols <- length(levels(factor(df.gnn.long$popGroup)))
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

ggplot(df.gnn.long, aes(x = fct_reorder(haplo, order), y = prob, fill = factor(popGroup))) +
    geom_col(aes(color = popGroup), size = 0.1) +
    facet_grid(~taxa, switch = "x", scales = "free", space = "free") +
    theme_minimal() + labs(x = "Individuals", y = "GNN") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    ggtitle(title) +
    theme(
        panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_text(angle=90),
        panel.grid = element_blank(),
        strip.text.x = element_text(angle=0),
        plot.title = element_text(hjust=0.5)
    ) + xlab(NULL) + scale_color_manual("Population", values=myColors) + scale_fill_manual("Population", values=myColors)

ggsave(snakemake@output[["png"]], width = 25, height = 6)
