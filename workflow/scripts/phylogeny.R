#!/usr/bin/env Rscript
# -*- mode: R -*-

suppressPackageStartupMessages(
    {
        require(tidyverse)
        require(ape)
        require(phyclust)
        require(ggfortify)
        require(ggtree)
    })

input <- '/home/hugo/projects/phage-evo-paper/results/single/pggb/distance_matrix_removed_ecoli_removed_phages.tsv'
output <- 'plots'

colfunc <- colorRampPalette(c("red", "yellow"))
phage.groups <- c('S2-55s', paste0('P',1:10))
phage.colors <- c('#faa0f4', colfunc(10))
names(phage.colors) <- phage.groups
  
y <- read.delim(input)
y %>% mutate(path.a=path.a, path.b=path.b, jaccard=jaccard, path.a.length=NULL, path.b.length=NULL, intersection=NULL, euclidean=NULL) %>% pivot_wider(names_from=path.b, values_from=jaccard) %>% replace(is.na(.), 0) -> y.dist
y.tree <- nj(as.dist(y.dist[, !names(y.dist) %in% c("path.a")]))

group.info <- split(y.tree$tip.label, gsub("#.+", "", y.tree$tip.label))
y.tree <- groupOTU(y.tree, group.info)
p <- ggtree(y.tree, branch.length = 'none', layout = 'daylight') + geom_tippoint(aes(color=group), size=2) + scale_color_manual('Passages' , values=phage.colors)
dev.off()

# ggsave(paste(output, "ggtree.passage.pdf", sep="."), height=40, width=9)
